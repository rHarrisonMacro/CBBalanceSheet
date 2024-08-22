function [ySim,OBCsAreActiveSim,diagnostics] = ...
    project_LSS_model_under_OC_subject_to_non_policy_constraints(...
    Model,y0,Shocks,OptPolInfo,nonPolOBCinfo,Options)
% This function projects model variables from an LSS MAPS model under
% optimal commitment, subject to constraints on both policy instruments and
% non-instrument variables.
%
% INPUTS
%   -> Model: MAPS LSS model
%   -> y0: 2*(nx-1)+1 vector of initial conditions for endogenous variables.
%   -> Shocks: structure containing fields:
%       - anticipated: nz*H matrix of anticipated shock values
%   -> OptPolInfo: structure containing standard information about the
%      optimal policy problem. It has fields:
%       - policyEqNames: single string or column cell string array
%       - policyShockMnems: single string or column cell string array (use
%         emmpty string ('') if there are no policy shocks to remove)
%       - instrumentMnems: single string or column cell string array
%       - instrumentWeights: column vector of loss weights for instruments
%       - objVarMnems: single string or column cell string array
%       - objVarWeights: column vector of loss weights for objective vars
%       - beta: discount factor
%   -> nonPolOBCinfo: structure containing information about the
%      constraints on non-instruments. These are in fields "active" and
%     "inactive", which both contain a field "VarMnems", containing the
%     mnemnomics of the variables that are constrained when constraints are
%     active or inactive. The field "active" also contains a "Constants"
%     field, which contains the values that the constrained variables must
%     satisfy when the constraints are active.
%   -> Options (optional): structure containing solution options for the
%      non-instrument constraints. Fields are:
%        - solutionUpdateMethod: string specifying method to be used to
%        update guesses of occasionally binding constraints.
%        - guessUpdateHandle: if solutionUpdateMethod is set to
%        "UserDefined", then a function handle must be supplied.
%
% OUPUTS
%   -> ySim: ny*H simulation of endogenous variables and Lagrange
%      multipliers on the structural model equations.
%   -> OBCsAreActiveSim: N*H matrix of logicals, indicating when
%      constraints are "active" during the simulation.
%   -> diagnostics: structure containing diagnostic information about the
%      non-instrument constraints in the simulation output. See 
%      "check_simulation_for_non_pol_OBCs.m" for details of the fields
%      within this structure.
%
% 
% This version: 15/12/2020
% Author(s): Richard Harrison

%% INTERROGATE INPUTS
if isfield(Shocks,'unanticipated')
    error(['The Shocks structure contains an "unanticipated" field.',...
        'This function produces a perfect foresight projection with ',...
        'anticipated shocks only.']);
end
if nargin<5 || nargin>6
    error('This function requires 5 or 6 inputs.');
end

%% SET OPTIONS
if nargin<6
    Options = struct;
end
Options = overlay_default_options(Options);

%% UNPACK OPTIONS
solutionUpdateMethod = Options.solutionUpdateMethod;
if strcmpi(solutionUpdateMethod,'userDefined')
    if ~isfield(Options,'guessUpdateHandle')
        error(['Function handle must be supplied within Options stucture',...
            ' if solution method "userDefined" is specified.']);
    end
    guessUpdateHandle = Options.guessUpdateHandle;
    solutionMethodIsUserDefined = true;
else
    solutionMethodIsUserDefined = false;
end
maxIter = Options.maxIter;
instBoundTol = Options.instBoundTol;
if isempty(Options.initialGuess)
    initialiseOBCsInactive = true;
else
    initialiseOBCsInactive = false;
end

%% EXTRACT SIMULATION AND SHOCK INFORMATION
H = size(Shocks.anticipated,2);

%% PRODUCE THE BASE INPUTS FOR OPTIMAL COMMITMENT POLICY
% There is some replication of the standard OC solution codes here,
% reflecting the fact that outputs of high-level functions are re-ordered
% to respect the original variable ordering in the model. This code needs
% to preserve the low-level information.

%% UNPACK MODEL
[HB,HC,HF,PSI,xMnems,zMnems,xEqNames] = unpack_model(...
    Model,{'HB','HC','HF','PSI','xMnems','zMnems','xEqNames'});

%% UNPACK THE OPTPOLINFO STRUCTURE
policyEqNames = OptPolInfo.policyEqNames;
policyShockMnems = OptPolInfo.policyShockMnems;
instrumentMnems = OptPolInfo.instrumentMnems;
objVarMnems = OptPolInfo.objVarMnems;
beta = OptPolInfo.beta;
objVarWeights = OptPolInfo.objVarWeights;
instrumentWeights = OptPolInfo.instrumentWeights;

%% PARTITION THE STRUCTURAL EQUATIONS
[HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztildeNonPolShk,~,instrumentInds,~,~,policyShockLogicals,~] = ...
    partition_structural_equations_for_optimal_discretion_solution(...
    HB,HC,HF,PSI,xMnems,zMnems,xEqNames,instrumentMnems,policyEqNames,...
    policyShockMnems,objVarMnems);

%% CREATE ENDOGENOUS VARIABLES AND INSTRUMENTS LOSS FUNCTION WEIGHTS
% These are needed to build the time-varying HCy matrices
[W,Q] = create_loss_function_weight_matrices(...
    xMnems,objVarMnems,objVarWeights,instrumentMnems,instrumentWeights);
OptPolInfo.W = W;
OptPolInfo.Q = Q;
OptPolInfo.instrumentInds = instrumentInds;

%% COMPUTE "BASELINE" OC MATRICES AND UPDATE OPTPOLINFO STRUCTURE
nz = size(PSI,2);
nxTilde = size(HtildeBxtilde,1);
PSItildeztilde = zeros(nxTilde,nz);
PSItildeztilde(:,~policyShockLogicals) = PSItildeztildeNonPolShk;
[By,~,~,HBy,~,HFy,PSIy,OptPolInfo,indicesInfo] = ...
    compute_OC_solution_using_Dennis_algorithm(...
    HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,beta,OptPolInfo);

%% EXTRACT SHOCK INFORMATION
nmu = size(OptPolInfo.Constraints.shadowShockMnems,1);
zMat = Shocks.anticipated;
zMat = [zMat; zeros(nmu,H)];
muInds = ((nz+1):(nz+nmu))';

%% EXPAND NON-POLICY OBC INFORMATION STRUCTURE
nonPolOBCinfo = ...
    add_model_info_to_non_pol_OBC_struct(nonPolOBCinfo,Model,OptPolInfo);
PSIdelta = nonPolOBCinfo.PSIdelta;

%% EXTRACT AND ANALYSE INFORMATION ABOUT NON-POLICY OBCS
N = size(nonPolOBCinfo.G,1);

%% INITIALISE INPUTS
if initialiseOBCsInactive
    OBCsAreActive = false(N,H);
else
    OBCsAreActive = Options.initialGuess;
    if size(OBCsAreActive,1)~=N || size(OBCsAreActive,2)~=H 
        error('Initial guess for OBC status must be an N*T matrix.')
    end
end
    
PSIhaty = [PSIy, [zeros(indicesInfo.nx,N); PSIdelta]];

%% EXTRACT INFORMATION ABOUT THE INSTRUMENT CONSTRAINTS
instInds = lookup_model_index_numbers(xMnems,...
    OptPolInfo.Constraints.instrumentMnems);
instCoeffs = repmat(OptPolInfo.Constraints.instrumentCoeffs,1,H);
instBounds = repmat(OptPolInfo.Constraints.constants,1,H);


%% GET READY FOR LOOP
converged = false;
iter = 1;
cyclingDetected = false;
activeOBCguesses = nan(N,H,maxIter);

%% EXECUTE THE LOOP
while iter<maxIter && ~converged && ~cyclingDetected  
    % Construct matrices for time-varying solution
    [HCmatrices,delta] = build_time_varying_model_matrices(nonPolOBCinfo.G,...
        nonPolOBCinfo.EqInds,OBCsAreActive,HtildeCxtilde,HtildeCr,W,Q,H,...
        indicesInfo);

    [ByMats,PHIhatymats,Fmats] = compute_time_varying_solution_matrices(...
        HCmatrices,HFy,HBy,PSIhaty,By);

    % Update shocks to be consistent with current OBC guess
	zhat = [zMat; delta];
        
    % Project using the time-varying matrices
    ySim = project_using_time_varying_solution(ByMats,PHIhatymats,Fmats,...
        zhat,y0);
    
    % Check policy OBCs and impose if necessary
    instPaths = ySim(instInds,:);
    instPathDiffs = instPaths.*instCoeffs - instBounds;
    instrumentConstraintsViolated = any(any(instPathDiffs<instBoundTol));
    if instrumentConstraintsViolated
        % Build matrices
        M = build_HP_matrices(ByMats,PHIhatymats,Fmats,muInds,instInds);
        % Invert
        mustar = invert_to_find_shadow_shocks_using_Holden_Paetz(...
            instPaths,OptPolInfo.Constraints.constants,M);
        muShks = reshape(mustar,nmu,H);
        % Adjust zhat shocks and project
        zhatWithMu = zhat;
        zhatWithMu(muInds,:) = muShks;
        ySim = project_using_time_varying_solution(ByMats,PHIhatymats,...
            Fmats,zhatWithMu,y0);
    end
    
    % Check validity of solution
    [simIsValid,OBCsAreActiveSim,diagnostics] = ...
        check_simulation_for_non_pol_OBCs(ySim,nonPolOBCinfo,OBCsAreActive);
  
    % Check convergence and update guess for OBCs
    if simIsValid
        disp(['Equilibrium found on iteration ' num2str(iter)]);
        converged = true;
    else
        if solutionMethodIsUserDefined
            OBCsAreActive = update_non_policy_OBC_guess(OBCsAreActive,...
                OBCsAreActiveSim,diagnostics,solutionUpdateMethod,...
                guessUpdateHandle);
        else
            OBCsAreActive = update_non_policy_OBC_guess(OBCsAreActive,...
                OBCsAreActiveSim,diagnostics,solutionUpdateMethod);
        end
    end
    
    activeOBCguesses(:,:,iter) = double(OBCsAreActive);
    
    %% CHECK FOR CYCLING
    if iter>1
        OBChistory = reshape(activeOBCguesses(:,:,1:iter),H*N,iter);
        previousTries = OBChistory(:,1:iter-1);
        latestGuess = OBChistory(:,iter);
        identicalGuesses = (previousTries==latestGuess);
        cyclingTest = all(identicalGuesses);
        [~,repeatGuessInds] = find(cyclingTest==true,1,'last');
        if ~isempty(repeatGuessInds)
            cyclingDetected = true;
        end
    end
    
    %% PREPARE FOR NEXT ITERATION
    iter = iter + 1;    
end

%% CHECK FOR CYCLING
if cyclingDetected && ~converged
    error(['Cycling detected on iteration ' num2str(iter-1)]);
end

%% CHECK CONVERGENCE
if ~converged
    error('Solution did not converge in maximum number of iterations.');
end

end


%% SUBFUNCTION TO OVERLAY DEFAULTS IF OPTIONS STRUCTURE DOES NOT CONTAIN THEM
function Options = overlay_default_options(Options)

if ~isfield(Options,'initialGuess')
    Options.initialGuess = [];
end

if ~isfield(Options,'solutionUpdateMethod')
    Options.solutionUpdateMethod = 'fullUpdateInactiveFirst';
end

if ~isfield(Options,'maxIter')
    Options.maxIter = 100;
end

if ~isfield(Options,'convTol')
    Options.convTol = 1e-06;
end

if ~isfield(Options,'instBoundTol')
    Options.instBoundTol = 1e-06;
end


end

%% SUBFUNCTION TO BUILD TIME-VARYING H^C_Y MATRIX
function [HCmatrices,delta] = ...
    build_time_varying_model_matrices(G,EqInds,OBCsAreActive,...
    HtildeCxtilde,HtildeCr,W,Q,T,indicesInfo)
%% INTERROGATE AND UNPACK INFORMATION ABOUT PROBLEM DIMENSIONS
nx = indicesInfo.nx;
ny = indicesInfo.ny;
nr = indicesInfo.nr;
rInds = indicesInfo.rInds;
lambdaInds = indicesInfo.lambdaInds;
xtildeInds = indicesInfo.xtildeInds;

HCmatrices = nan(ny,ny,T);
HCyt = zeros(ny,ny);

%% BUILD H^C_{Y} MATRICES
for t=1:T
    
    HtildeCxtildet = HtildeCxtilde;
    if any(OBCsAreActive(:,t))
        eqsToReplace = EqInds(OBCsAreActive(:,t));
        HtildeCxtildet(eqsToReplace,:) = G(OBCsAreActive(:,t),:);
    end

    % Load first order conditions for instrument:
    HCyt(1:nr,rInds) = Q;
    HCyt(1:nr,lambdaInds) = -HtildeCr';
    % Load first order conditions for non-policy variables:
    HCyt(nr+1:nx,xtildeInds) = W;
    HCyt(nr+1:nx,lambdaInds) = -HtildeCxtildet';
    % Load first order conditions with respect to multipliers:
    HCyt(nx+1:ny,xtildeInds) = HtildeCxtildet;
    HCyt(nx+1:ny,rInds) = HtildeCr;
    
    HCmatrices(:,:,t) = HCyt;
    
end

%% BUILD DUMMY SHOCK VALUES
delta = double(OBCsAreActive);

end

%% SUBFUNCTION TO ITERATE BACKWARDS AND SOLVE FOR TIME-VARYING SOLUTION MATRICES
function [ByMats,PHIhatymats,Fmats] = compute_time_varying_solution_matrices(...
        HCmats,HFy,HBy,PSIhaty,By)

[ny,~,T] = size(HCmats);    
nz = size(PSIhaty,2);    

ByMats = nan(ny,ny,T);    
eyeny = eye(ny);
PHIhatymats = nan(ny,nz,T);   
Fmats = nan(ny,ny,T,T);

Bytplus1 = By;
for t=T:-1:1
    invMat = eyeny/(HFy*Bytplus1 + HCmats(:,:,t));
    Byt = -invMat*HBy;
    ByMats(:,:,t) = Byt;
    PHIhatymats(:,:,t) = invMat*PSIhaty;
    % F matrices
    Upsilont = -invMat*HFy;
    for h = 1:T-t
        if h==1
            Fmats(:,:,t,h) = Upsilont;
        else
            Fmats(:,:,t,h) = Upsilont*Fmats(:,:,t+1,h-1);
        end
    end
    % Prepare for next iteration
    Bytplus1 = Byt;
end

end

%% SUB-FUNCTION TO PROJECT WITH TIME-VARYING LAWS OF MOTION
function     ySim = project_using_time_varying_solution(ByMats,PHIhatymats,Fmats,...
        zhat,y0)

    yLag = y0;
    T = size(zhat,2);
    ySim = nan(size(ByMats,1),T);
    for t=1:T
        shockEffects = PHIhatymats(:,:,t)*zhat(:,t);
        for h = 1:T-t
            shockEffects = shockEffects + ...
                Fmats(:,:,t,h)*PHIhatymats(:,:,t+h)*zhat(:,t+h);
        end
        yt = ByMats(:,:,t)*yLag + shockEffects;
        ySim(:,t) = yt;
        yLag = yt;
    end
    
end

%% SUB-FUNCTION TO BUILD HOLDEN-PAETZ MATRICES 
function M = build_HP_matrices(ByMats,PHIhatymats,Fmats,muInds,instInds)

%% DEFINE SELECTOR MATRIX TO ISOLATE IMPACTS ON INSTRUMENTS
ny = size(ByMats,1);
nmu = size(muInds,1);
Stau = zeros(nmu,ny);
for imu = 1:nmu
    Stau(imu,instInds(imu)) = 1;
end

%% INITIALISE MATRIX
T =  size(ByMats,3);
M = nan(nmu*T,nmu*T);
Rlag = zeros(ny,nmu*T);

%% BUILD MATRIX
for t = 1:T
    omega = zeros(ny,nmu*T);
    PHImut = PHIhatymats(:,muInds,t);
    omega(:,nmu*(t-1)+1:nmu*t) = PHImut;
    for h=1:T-t
        PHImuth = PHIhatymats(:,muInds,t+h);
        omega(:,nmu*(t+h-1)+1:nmu*(t+h)) = Fmats(:,:,t,h)*PHImuth;
    end
    Rt = ByMats(:,:,t)*Rlag + omega;
    M(nmu*(t-1)+1:nmu*t,:) = Stau*Rt;
    Rlag = Rt;
end

end