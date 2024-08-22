function [xSim,OBCsAreActiveSim,diagnostics] = ...
    project_LSS_model_under_OD_subject_to_non_policy_constraints(...
    Model,x0,Shocks,OptPolInfo,nonPolOBCinfo,Options)
% This function projects model variables from an LSS MAPS model under
% optimal discretion, subject to constraints on both policy instruments and
% non-instrument variables.
%
% INPUTS
%   -> Model: MAPS LSS model
%   -> x0: nx*1 vector of initial conditions for endogenous variables.
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
% NOTES
%   -> This function can be viewed as a "wrapper" for the lower level
%      functions used to compute discretionary equilibria subject to
%      instrument bounds. It could be usefully (and easily) extended to 
%      include additional functionality allowing the user to pass more 
%      solution options to the lower level functions.
% 
% This version: 15/12/2020
% Author(s): Richard Harrison.

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
maxIter = Options.maxIter;
reportProgress = Options.reportProgress;

if isempty(Options.initialGuess)
    initialiseOBCsInactive = true;
else
    initialiseOBCsInactive = false;
end

%% EXTRACT SIMULATION AND SHOCK INFORMATION
H = size(Shocks.anticipated,2);

%% SOLVE MODEL UNDER OPTIMAL DISCRETIONARY POLICY
% Note that additional functionality could be added to allow this function
% to pass solution options to the OD solver called here. At present, the
% solution is carried out with the default options.
[BOD,~,HBOD,HCOD,HFOD,PSIOD,~,ODinfo,Vxtildextilde] = ...
    solve_LSS_model_under_optimal_discretion(Model,OptPolInfo); 

%% CONTROL OD OPTIONS WITH HARD-WIRED LOGICAL
% As noted above, functionality could be added to facilitate a richer set 
% of options being passed to the function. This logical determines whether
% the indicator guess for the *policy* OBCs is initialised to zero or a
% user-supplied function handle.
initialiseWithZeros = true;

%% UNPACK MNEMONIC, LOGICALS AND INDEX INFORMATION FROM ODINFO STRUCTURE
instrumentMnems = ODinfo.instrumentMnems;
instrumentLogicals = ODinfo.instrumentLogicals;
instrumentInds = ODinfo.instrumentInds;
policyEqLogicals = ODinfo.policyEqLogicals;
policyShockLogicals = ODinfo.policyShockLogicals;

%% UNPACK LOSS FUNCTION INFO
beta = ODinfo.beta;
W = ODinfo.W;
Q = ODinfo.Q;

%% PARTITION B & x0
Bxtildextilde = BOD(~instrumentLogicals,~instrumentLogicals);
Brxtilde = BOD(instrumentInds,~instrumentLogicals);
xtilde0 = x0(~instrumentLogicals);

%% PARTITION STRUCTURAL EQUATION MATRICES
HtildeBxtilde = HBOD(~policyEqLogicals,~instrumentLogicals);
HtildeCxtilde = HCOD(~policyEqLogicals,~instrumentLogicals);
HtildeCr = HCOD(~policyEqLogicals,instrumentInds);
HtildeFxtilde = HFOD(~policyEqLogicals,~instrumentLogicals);
HtildeFr = HFOD(~policyEqLogicals,instrumentInds);
PSItildeztilde = PSIOD(~policyEqLogicals,~policyShockLogicals);

%% UNPACK AND TRANSLATE CONSTRAINT INFO
% Constraints are of the form S*r >= b
[S,b] = create_instrument_bound_constraint_matrices(...
    instrumentMnems,ODinfo.Constraints);

%% COMPUTE DIM OF ENDOGENOUS VARIABLES, INSTRUMENTS, CONSTRAINTS & SHOCKS
nx = size(x0,1);
nmu = size(b,1);

%% EXTRACT SHOCK INFORMATION
zMat = Shocks.anticipated;
zMat = zMat(~policyShockLogicals,:);

%% EXPAND NON-POLICY OBC INFORMATION STRUCTURE
nonPolOBCinfo = ...
    add_model_info_to_non_pol_OBC_struct(nonPolOBCinfo,Model,OptPolInfo);
PSIdelta = nonPolOBCinfo.PSIdelta;

%% EXTRACT AND ANALYSE INFORMATION ABOUT NON-POLICY OBCS
N = size(nonPolOBCinfo.G,1);

%% BUILD SHOCK LOADING INFORMATION
PSIhat = [PSItildeztilde, PSIdelta];

%% INITIALISE INPUTS
if initialiseOBCsInactive
    OBCsAreActive = false(N,H);
else
    OBCsAreActive = Options.initialGuess;
    if size(OBCsAreActive,1)~=N || size(OBCsAreActive,2)~=H 
        error('Initial guess for OBC status must be an N*T matrix.')
    end
end

%% GET READY FOR LOOP
converged = false;
iter = 1;
xSim = nan(nx,H);
cyclingDetected = false;
activeOBCguesses = nan(N,H,maxIter);
% TO RECORD ALL ITERATION ATTEMPTS (FOR DEBUGGING)
allSims = cell(maxIter,1);

%% EXECUTE THE LOOP
while iter<=maxIter && ~converged && ~cyclingDetected
    % Construct matrices for time-varying solution
    [HCmatrices,delta] = build_time_varying_model_matrices(nonPolOBCinfo.G,...
        nonPolOBCinfo.EqInds,OBCsAreActive,HtildeCxtilde);

    % Update shocks to be consistent with current OBC guess
	zhat = [zMat; delta];
    
    if initialiseWithZeros
        J0 = zeros(nmu,nmu,H);
    else
        try                                                                 %#ok<UNRCH>
            J0 = initialisationFuncHandle(xtilde0,zhat,H);
        catch InitialisationE
            warning(['Function to initialise indicator for BPY ',...
                'failed with the following message ''',...
                InitialisationE.message, ''' ... initialising ',...
                'with zeros instead']);
            J0 = zeros(nmu,nmu,H);
        end
    end
        
    % Find constrained OD solution conditional on guess for non-policy OBCs
    [xtildeSim,rSim] = ...
            compute_optimal_discretion_transition_using_BPY(xtilde0,...
            zhat,HtildeBxtilde,HCmatrices,HtildeCr,...
            HtildeFxtilde,HtildeFr,PSIhat,beta,W,Q,S,b,J0,...
            Vxtildextilde,Bxtildextilde,Brxtilde); %,BPYoptions);  
	xSim(~instrumentLogicals,:) = xtildeSim;
    xSim(instrumentLogicals,:) = rSim;
    
    % Check validity of solution
    [simIsValid,OBCsAreActiveSim,diagnostics] = ...
        check_simulation_for_non_pol_OBCs(xSim,nonPolOBCinfo,OBCsAreActive);
  
    % Record simulation & assumptions under which it was made (debugging)
    allSims{iter} = xSim;
    activeOBCguesses(:,:,iter) = double(OBCsAreActive);
    
    % Check convergence and update guess for OBCs
    if simIsValid
        disp(['Equilibrium found on iteration ' num2str(iter)]);
        converged = true;
    else
        switch solutionUpdateMethod
            case 'userDefined'
                if ~isfield(Options,'guessUpdateHandle')
                    error(['Function handle to update OBC indicators ',...
                        'must be provided if the "userDefined" method ',...
                        'is selected.']);
                end
                OBCsAreActive = update_non_policy_OBC_guess(OBCsAreActive,...
                    OBCsAreActiveSim,diagnostics,solutionUpdateMethod,...
                    Options.guessUpdateHandle);
            otherwise
                OBCsAreActive = update_non_policy_OBC_guess(OBCsAreActive,...
                    OBCsAreActiveSim,diagnostics,solutionUpdateMethod);
        end
    end

    % Check for cycling
    if iter>2
        OBChistory = reshape(activeOBCguesses(:,:,1:iter),H*N,iter);
        cyclingDetected = check_OBC_guesses_for_cycling(OBChistory);
    end
    
    if reportProgress
        disp(['Iteration ' num2str(iter) ' complete.']);
    end
    
    % Prepare for next iteration
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

if ~isfield(Options,'reportProgress')
    Options.reportProgress = false;
end

end

%% SUBFUNCTION TO BUILD TIME-VARYING H^C_Y MATRIX
function [HCmatrices,delta] = ...
    build_time_varying_model_matrices(G,EqInds,OBCsAreActive,...
    HtildeCxtilde)
%% INTERROGATE AND UNPACK INFORMATION ABOUT PROBLEM DIMENSIONS
nxtilde = size(HtildeCxtilde,1);
H = size(OBCsAreActive,2);

%% INITIALISE OUTPUT
HCmatrices = nan(nxtilde,nxtilde,H);

%% BUILD H^C_{\tilde x,t} MATRICES
for t=1:H
    HtildeCxtildet = HtildeCxtilde;
    if any(OBCsAreActive(:,t))
        eqsToReplace = EqInds(OBCsAreActive(:,t));
        HtildeCxtildet(eqsToReplace,:) = G(OBCsAreActive(:,t),:);
    end
    HCmatrices(:,:,t) = HtildeCxtildet;
end

%% BUILD DUMMY SHOCK VALUES
delta = double(OBCsAreActive);

end
