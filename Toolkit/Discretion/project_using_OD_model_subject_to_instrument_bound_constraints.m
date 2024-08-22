function [xf,Jf,muf,xfe,Jfe,mufe,S,b] = ...
    project_using_OD_model_subject_to_instrument_bound_constraints(...
    ODinfo,HBOD,HCOD,HFOD,PSIOD,BOD,Vxtildextilde,x0,Shocks,Htail,Options)
% Projects model solved under discretion subject to instrument constraints.
%
% INPUTS:
%   -> ODinfo: structure of info about OD problem and constraints:
%       - policyEqNames: single string or column cell string array
%       - policyShockMnems: single string or column cell string array (use
%         emmpty string ('') if there are no policy shocks to remove)
%       - instrumentMnems: single string or column cell string array
%       - instrumentWeights: column vector of loss weights for instruments
%       - objVarMnems: single string or column cell string array
%       - objVarWeights: column vector of loss weights for objective vars
%       - beta: discount factor
%       - policyEqLogicals: vector of logicals that locate policy eqs
%       - policyEqInds: vector of indices that locate policy eqs
%       - policyShockLogicals: vector of logicals locating policy shocks
%       - policyShockInds: vector of indices locating policy shocks
%       - instrumentLogicals: vector of logicals locating the instruments
%       - instrumentInds: vector of indices that locate the instruments
%       - Q: instruments weights matrix
%       - W: target weights matrix
%       - Constraints: (optional/application specific) sructure describing
%         S*r >= b:
%           - instrumentMnems: single string or column cell string array
%           - instrumentCoeffs: numeric column vector of coeffs on 
%           - constants: numeric column vector of constraint constants
%   -> HBOD: loadings on lagged variables in structural eqs
%   -> HCOD: loadings on contemporaneous variables in structural eqs
%   -> HFOD: loadings on one-period ahead variables in structural eqs
%   -> PSIOD: loadings on shocks in structural eqs
%   -> BOD: nx*nx transtion matrix from unconstrained Dennis solution
%   -> Vxtildextilde: nxtilde*nxtilde quadratic coefficients on 
%      non-instrument endogenous variables in period ahead loss function 
%      from unconstrained (Dennis) solution
%   -> x0: nx*1 vector of initial conditions for model variables
%   -> Shocks: struture with: 
%       - anticipated (optional): nz*H matrix of shocks
%       - unanticipated (optional): nz*H matrix of shocks
%   -> Htail: additional number of periods over which to compute 
%      expectations
%   -> Options (optional): 
%       - initialisationFuncHandle (optional): function to use to 
%         initialise constraint binding indicator function for BPY
%       - BPY (optional): options for the BPY algorithm (if no input, then
%         defaults from the function called are used)
%           - tol: tolerance for convergence
%           - maxIter: maximum number before giving up
%           - reportProgress: false/true
%           - useBestOnNonConvergence: false/true
%
% OUTPUTS
%   -> xf: nx*H matrix of projected model variables
%   -> Jf: nmu*nmu*H matrix of indicators for when each of the constraints
%      is binding
%   -> muf: nmu*H matrix of projected Lagrange multipliers
%   -> xfe: nx*(H+Htail)*H matrix of endog var expectations
%   -> Jfe: nmu*nmu*(H+Htail)*H matrix of indicators for when each of the 
%      constraints is binding in expectation
%   -> mufe: nmu*(H+Htail)*H matrix of shadow shock expectations
%   -> S: nmu*nr matrix of coefficients in the bound constraints
%   -> b: nmu*1 vector of constants in the bound constraints
%
% DETAILS:
%   -> This function is designed to support three main use cases: a) an 
%      ODPP with anticipated shocks as judgements; b) impulse responses and
%      other projections which could be a mix of anticipated and 
%      unanticpated shocks; c) stochastic simulation, which would be a 
%      sequence of unanticipated shocks. The function is designed to do the
%      minimum required, which depends on which of the two sets off shocks
%      is present.
%   -> The initialisationFuncHandle must accept the initial condition for
%      the non-instrumente endogenous variables and a set of non-policy
%      anticipated shocks as the first two inputs.  All other inputs 
%      necessary for construction of the guess must be passed to the
%      function anonymously at the time of the creation of the handle.
%
% NOTES:   
%   -> The steady state regime that is returned to over the transition is 
%      one in which none of the instrument constraints is binding.
%   -> The BPY algorithm is initialised with a guess that none of the
%      instrument constraints is binding.  It would be possible to extend
%      this function to permit alternative initialisations (e.g. based on
%      some pre-programmed alternatives implemented in sub-functions or a 
%      general function handle implementation).
%   -> Note that error handling is kept to a minimum.
%
% This version: 22/03/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 10
    error('This function requires 10 or 11 inputs')
end

%% HANDLE OPTIONAL OPTIONS INPUT
initialiseWithZeros = true;
BPYoptions = struct;
if nargin > 10
    check_field_names_in_structure(...
        Options,'optional options input','',...
        {'initialisationFuncHandle';'BPY'});
    if isfield(Options,'initialisationFuncHandle')
        initialiseWithZeros = false;
        initialisationFuncHandle = Options.initialisationFuncHandle;
    end
    if isfield(Options,'BPY')
        BPYoptions = Options.BPY;
    end        
end

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

%% DETERMINE NATURE OF PROJECTION & COMPUTE FORECAST HORIZON FROM SHOCKS
projectionContainsAnticipatedShocks = false;
projectionContainsUnanticipatedShocks = false;
if isfield(Shocks,'anticipated')
    projectionContainsAnticipatedShocks = true;
    H = size(Shocks.anticipated,2);
end
if isfield(Shocks,'unanticipated')
    projectionContainsUnanticipatedShocks = true;
    H = size(Shocks.unanticipated,2);
end
if ~projectionContainsAnticipatedShocks && ...
        ~projectionContainsUnanticipatedShocks
    error(['This function requires a set of anticipated and/or ',...
        'unanticipated shocks (i.e. can''t be neither)']);
end

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
nxtilde = size(W,1);
nr = size(Q,1);
nmu = size(b,1);
nz = size(PSIOD,2);

%% CONSTRUCT RECURSIVE SHOCKS
% By the nature of the problem, it is necessary to recompute expectations
% (i.e. for when the constraints will bind in future) whenever agents' 
% information sets change (i.e. whenever there is an unanticipated shock).
% To simplify below we construct a recursive set of anticipated shocks, 
% zfe, where zfe(:,:,h) is the set of shocks applying in expectation in 
% period h. zfe(:,s,h) references expected shocks in period s, where 
% zfe(:,s,h)=NaN for all s<h and zfe(:,s,h)=0 for all s>H. The number of 
% periods over which it is necessary to compute expectations depends on 
% whether or not the problem contains anticipated shocks and on the user 
% choice in the input Htail. If the projection does contain anticipated 
% shocks then, in each period, it is necessary to compute expectations all 
% the way out to H+Htail. If it does not, then we only need to compute 
% expectations out to 1+Htail.
zfe = NaN*ones(nz,H+Htail,H);
for h = 1:H
    zfe(:,h:h+Htail,h) = zeros(nz,Htail+1);
end
if projectionContainsAnticipatedShocks
    for h = 1:H
        zfe(:,h:H+Htail,h) = zeros(nz,H+Htail-h+1);
        zfe(:,h:H,h) = Shocks.anticipated(:,h:H);
    end
end
if projectionContainsUnanticipatedShocks
    for h = 1:H
        if any(Shocks.unanticipated(:,h)~=0)
            zfe(:,h,h) = zfe(:,h,h)+Shocks.unanticipated(:,h);
        end
    end
end

%% CREATE A RECRUSIVE SET OF SHOCKS WITHOUT THE POLICY SHOCKS
ztildefe = zfe(~policyShockLogicals,:,:);

%% COMPUTE MAXIMAL HORIZON FOR EXPECTATIONS
if projectionContainsAnticipatedShocks
    Hfe = H+Htail;
else
    Hfe = 1+Htail;
end

%% INITIALISE RECUESIVE EXPECTATIONS AND PROJECTION RESULT MATRICES
xfe = NaN*ones(nx,H+Htail,H);
xtildefe = NaN*ones(nxtilde,H+Htail,H);
rfe = NaN*ones(nr,H+Htail,H);
mufe = NaN*ones(nmu,H+Htail,H);
Jfe = NaN*ones(nmu,nmu,H+Htail,H);
xf = NaN*ones(nx,H);
xtildef = NaN*ones(nxtilde,H);
rf = NaN*ones(nr,H);
muf = NaN*ones(nmu,H);
Jf = NaN*ones(nmu,nmu,H);

%% MAIN ALGORITHM
% The algorithm is designed to compute the mimimum necessary.  It always
% computes a projection for h=1, but then only updates that projection if
% new information (unanticipated shocks) arrive.  Note also that, if the
% projection contains anticipated shocks, then the horizon over which 
% expected shadow shocks is computed is reduced as we move through the 
% projection horizon.
xtilde0h = xtilde0;
Hfeh = Hfe;
for h = 1:H
    %% PROJECT IN PERIOD 1 & IN ALL PERIODS IF UNANTICIPATED SHOCKS PRESENT
    if h==1 || projectionContainsUnanticipatedShocks
        %% EXTRACT SHOCKS TO APPLY
        ztildefeh = ztildefe(:,h:h+Hfeh-1,h);
        %% INITIALISE BINDING CONSTRAINT INDICATOR
        if initialiseWithZeros
            J0 = zeros(nmu,nmu,Hfeh);
        else
            try
                J0 = initialisationFuncHandle(xtilde0h,ztildefeh,Hfeh);
            catch InitialisationE
                warning(['Function to initialise indicator for BPY ',...
                    'failed with the following message ''',...
                    InitialisationE.message, ''' ... initialising ',...
                    'with zeros instead']);
                J0 = zeros(nmu,nmu,Hfeh);
            end
        end
        %% COMPUTE A PROJECTION
        [xtildefe(:,h:h+Hfeh-1,h),rfe(:,h:h+Hfeh-1,h),...
            mufe(:,h:h+Hfeh-1,h),Jfe(:,:,h:h+Hfeh-1,h)] = ...
            compute_optimal_discretion_transition_using_BPY(xtilde0h,...
            ztildefeh,HtildeBxtilde,HtildeCxtilde,HtildeCr,...
            HtildeFxtilde,HtildeFr,PSItildeztilde,beta,W,Q,S,b,J0,...
            Vxtildextilde,Bxtildextilde,Brxtilde,BPYoptions);       
        %% STORE THE RESULT
        xtildef(:,h:H) = xtildefe(:,h:H,h);        
        rf(:,h:H) = rfe(:,h:H,h);
        muf(:,h:H) = mufe(:,h:H,h);
        Jf(:,:,h:H) = Jfe(:,:,h:H,h); 
    end
    %% UPDATE INITIAL CONDITION
    xtilde0h = xtildef(:,h);
    %% ADJUSTED HORIZON FOR FORECASTS IN ANTICIPATED SHOCKS CASE
    if projectionContainsAnticipatedShocks
        Hfeh = Hfeh-1;
    end
end

%% CONSTRUCT FULL ENDOGENOUS VARIABLE PROJECTION MATRICES
xfe(~instrumentLogicals,:,:) = xtildefe;
xfe(instrumentInds,:,:) = rfe;
xf(~instrumentLogicals,:) = xtildef;
xf(instrumentInds,:) = rf;

end