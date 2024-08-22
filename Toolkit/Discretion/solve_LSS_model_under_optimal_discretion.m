function [BOD,PHIOD,HBOD,HCOD,HFOD,PSIOD,FODhandle,ODinfo,Vxtildextilde,...
    Theta,zeta] = solve_LSS_model_under_optimal_discretion(...
    Model,ODinfo,Options) 
% Solves MAPS LSS model under OD using a MAPS implementation of Dennis.
% The solution is of the same form as in standard MAPS models with the
% exception that the impact of anticipated shocks cannot be summarised 
% using a single matrix.
%
% INPUTS:   
%   -> Model: MAPS LSS model
%   -> ODinfo: structure of info required to solve the OD problem
%       - policyEqNames: single string or column cell string array
%       - policyShockMnems: single string or column cell string array (use
%         emmpty string ('') if there are no policy shocks to remove)
%       - instrumentMnems: single string or column cell string array
%       - instrumentWeights: column vector of loss weights for instruments
%       - objVarMnems: single string or column cell string array
%       - objVarWeights: column vector of loss weights for objective vars
%       - beta: discount factor
%   -> Options (optional): structure
%       - algorithmIsDennis (optional): true/false
%       - Initialisation (optional): Initialisation choice for B/V
%           - source (optional): 'setToZeros', 'useInputModel', 
%             'setToUserInput'
%           - Binit: nx*nx matrix for `setToUserInput' option
%           - Vinit: nxtilde*nxtilde matrix for 
%             `setToUserInput' option if algorithmIsDennis = false
%       - Algorithm (optional): structure
%           - tol (optional): tolerance for convergence
%           - maxIter (optional): maximum number of iterations
%
% OUTPUTS:  
%   -> BOD: state transition for model solved under optimal discretion
%   -> PHIOD: loadings on shocks in optimal discretion solution
%   -> HBOD: loadings on lagged variables in structural eqs
%   -> HCOD: loadings on contemporaneous variables in structural eqs
%   -> HFOD: loadings on one-period ahead variables in structural eqs
%   -> PSIOD: loadings on shocks in structural eqs
%   -> FODhandle: function handle that can be used to compute the
%      horizon-varying anticipated shock loadings for a user choice of
%      horizon (and treatment of the intervening periods and starting point 
%      -- i.e. does counting begin at F_1 or F_0) 
%   -> ODinfo: input structure updated with:
%       - policyEqLogicals: vector of logicals that locate policy eqs
%       - policyEqInds: vector of indices that locate policy eqs
%       - policyShockLogicals: vector of logicals locating policy shocks
%       - policyShockInds: vector of indices locating policy shocks
%       - instrumentLogicals: vector of logicals locating the instruments
%       - instrumentInds: vector of indices that locate the instruments
%       - Q: instruments weights matrix
%       - W: target weights matrix
%   -> Vxtildextilde: coefficients on quadratic of contemporaneous 
%      endogenous variables in period-ahead losses
%   -> Theta: "solved out" coefficient on contemporaneous endogenous vars 
%      in structural equations 
%      (HtildeCx+HtildeFxtilde*Bxtildextilde+HtildeFr*Brxtilde)
%   -> zeta: coefficient common to loadings on all variables in the
%      targeting rule with the constraint substituted in
%
% DETAILS: 
%   -> This function solves for optimal discretionary policy in the LQ
%      framework.
%   -> It compiles the "full" set of model matrices, including the solution
%      and structural equations.
%   -> It also provides a function handle that can be used to compute a
%      time-varying forward loading matrix for use in projections.
%   -> There are two alternative algorithms that can be used: the Dennis
%      algorithm (which is based on policy function or law of motion 
%      iteration); a value function iteration algorithm, which should
%      produce equivalent results but has different numerical properties.
%
% NOTES: 
%   -> Error handling is kept to a minimum, though note that there is some
%      checking of the compatibility of the inputs with each other and with
%      the underlying assumptions of the algorithms in the two input
%      preparation functions.
%
% This version: 06/03/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 2
    error('This function requires 2 or 3 inputs')
end

%% HANDLE OPTIONAL OPTIONS INPUT
if nargin < 3
    Options = struct;
end

%% UNPACK AND VALIDATE OPTIONS INPUT
% If no options passed in as input or a field is missing, defaults (set in
% the subfunction) are: algorithmIsDennis = true; 
% initialisationSource = 'useInputModel' for algorithmIsDennis = true and 
% initialisationSource = 'setToZeros' for algorithmIsDennis = false; 
% AlgorithmOpts = struct; (which will mean that default algorithm options 
% are used in the specific algorithm called below.
[algorithmIsDennis,initialisationSource,InitialisationData,...
    AlgorithmOpts] = unpack_and_validate_solution_options_structure(...
    Options);

%% UNPACK MODEL
[HB,HC,HF,PSI,B,xMnems,zMnems,xEqNames] = unpack_model(...
    Model,{'HB','HC','HF','PSI','B','xMnems','zMnems','xEqNames'});

%% UNPACK THE INFO STRUCTURE
policyEqNames = ODinfo.policyEqNames;
policyShockMnems = ODinfo.policyShockMnems;
instrumentMnems = ODinfo.instrumentMnems;
instrumentWeights = ODinfo.instrumentWeights;
objVarMnems = ODinfo.objVarMnems;
objVarWeights = ODinfo.objVarWeights;
beta = ODinfo.beta;

%% PARTITION THE STRUCTURAL EQUATIONS
[HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,instrumentLogicals,instrumentInds,policyEqLogicals,...
    policyEqInds,policyShockLogicals,policyShockInds] = ...
    partition_structural_equations_for_optimal_discretion_solution(...
    HB,HC,HF,PSI,xMnems,zMnems,xEqNames,instrumentMnems,policyEqNames,...
    policyShockMnems,objVarMnems);

%% CREATE ENDOGENOUS VARIABLES AND INSTRUMENTS LOSS FUNCTION WEIGHTS
[W,Q] = create_loss_function_weight_matrices(...
    xMnems,objVarMnems,objVarWeights,instrumentMnems,instrumentWeights);

%% UPDATE THE ODINFO STRUCTURE WITH THE LOGICALS AND LOSS FUNCTION WEIGHTS
ODinfo.instrumentLogicals = instrumentLogicals;
ODinfo.instrumentInds = instrumentInds;
ODinfo.policyEqLogicals = policyEqLogicals;
ODinfo.policyEqInds = policyEqInds;
ODinfo.policyShockLogicals = policyShockLogicals;
ODinfo.policyShockInds = policyShockInds;
ODinfo.Q = Q;
ODinfo.W = W;

%% COMPUTE DIMENSIONS OF OD MODEL FOR USE BELOW
nx = size(B,1);
nz = size(PSI,2);
nr = size(Q,1);
nxtilde = size(W,1);

%% INITIALISE SOLUTION
switch initialisationSource
    case 'setToZeros'
        B0 = zeros(nx,nx);
        if ~algorithmIsDennis
            V0xtildextilde = zeros(nxtilde,nxtilde);
        end
    case 'useInputModel'
        B0 = B;
        if ~algorithmIsDennis
            % NB: this is a slightly odd initialisation
            B0xtildextilde = B0(~instrumentLogicals,~instrumentLogicals);
            B0rxtilde = B0(instrumentInds,~instrumentLogicals);
            V0xtildextilde = ...
                compute_infinite_sum_using_doubling_algorithm(...
                beta*B0xtildextilde',...
                B0xtildextilde'*W*B0xtildextilde+B0rxtilde'*Q*B0rxtilde,...
                B0xtildextilde);
        end
    case 'setToUserInput'
        [B0,V0xtildextilde] = ...
            validate_and_unpack_user_input_initialisation(...
            InitialisationData,algorithmIsDennis,nx,nxtilde);
end
B0xtildextilde = B0(~instrumentLogicals,~instrumentLogicals);
B0rxtilde = B0(instrumentInds,~instrumentLogicals);

%% CALL FUNCTION TO RUN DENNIS ALGORITHM
if algorithmIsDennis
    [BODxtildextilde,BODrxtilde,PHIODxtildeztilde,PHIODrztilde,...
        HhatCxtilde,HhatCr,Vxtildextilde,Theta,zeta,Deltar] = ...
        compute_optimal_discretion_solution_using_Dennis_algorithm(...
        HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
        PSItildeztilde,beta,W,Q,B0xtildextilde,B0rxtilde,AlgorithmOpts);
end

%% CALL FUNCTION TO RUN VALUE FUNCTION ITERATION ALGORITHM
if ~algorithmIsDennis
    [BODxtildextilde,BODrxtilde,PHIODxtildeztilde,PHIODrztilde,...
        HhatCxtilde,HhatCr,Vxtildextilde,Theta,zeta,Deltar] = ...
        compute_optimal_discretion_solution_using_VFI(...
        HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
        PSItildeztilde,beta,W,Q,V0xtildextilde,B0xtildextilde,...
        B0rxtilde,AlgorithmOpts);
end

%% COMPILE THE SOLUTION MATRICES
BOD = zeros(nx,nx);
PHIOD = zeros(nx,nz);
BOD(~instrumentLogicals,~instrumentLogicals) = BODxtildextilde;
BOD(instrumentInds,~instrumentLogicals) = BODrxtilde;
PHIOD(~instrumentLogicals,~policyShockLogicals) = PHIODxtildeztilde;
PHIOD(instrumentInds,~policyShockLogicals) = PHIODrztilde;

%% REPLACE THE ORIGINAL INSTRUMENT EQUATIONS IN THE STRUCTURAL EQUATIONS
HBOD = zeros(nx,nx);
HCOD = zeros(nx,nx);
HFOD = zeros(nx,nx);
PSIOD = zeros(nx,nz);
HBOD(~policyEqLogicals,~instrumentLogicals) = HtildeBxtilde;
HCOD(~policyEqLogicals,~instrumentLogicals) = HtildeCxtilde;
HCOD(~policyEqLogicals,instrumentInds) = HtildeCr;
HCOD(policyEqInds,~instrumentLogicals) = HhatCxtilde;
HCOD(policyEqInds,instrumentInds) = HhatCr;
HFOD(~policyEqLogicals,~instrumentLogicals) = HtildeFxtilde;
HFOD(~policyEqLogicals,instrumentInds) = HtildeFr;
PSIOD(~policyEqLogicals,~policyShockLogicals) = PSItildeztilde;

%% CREATE THE LOSS FUNCTION WEIGHTING MATRIX FOR ANTICIPATED SHOCKS FORMULA
Omega = zeros(nx,nx);
Omega(instrumentInds,instrumentInds) = Q;
Omega(~instrumentLogicals,~instrumentLogicals) = W+beta*Vxtildextilde;

%% CREATE THE PRIVATE SECTOR & POLICYMAKER FORWARD LOADINGS MATRICES
% Denoted as \mathcal{F}^{ps} and \mathcal{F}^{pol} in the derivation.
Fps = zeros(nx,nx);
Fpol = zeros(nx,nx);
% Define inverse matrices for repeated usage
ThetaInv = eye(nxtilde,nxtilde)/Theta;
DeltarInv = eye(nr,nr)/Deltar;
% Define \Delta_{F} matrices
DeltaFpsrxtilde = -zeta*HtildeFxtilde;
DeltaFpsrr = -zeta*HtildeFr;
DeltaFpolrxtilde = (ThetaInv*HtildeCr)'*beta*BODxtildextilde';
DeltaFpolrr = (ThetaInv*HtildeCr)'*beta*BODrxtilde';
% Define F^{ps} matrices
Fps(~instrumentLogicals,~instrumentLogicals) = ...
    -ThetaInv*(HtildeFxtilde+HtildeCr*DeltarInv*DeltaFpsrxtilde);
Fps(~instrumentLogicals,instrumentInds) = ...
    -ThetaInv*(HtildeFr+HtildeCr*DeltarInv*DeltaFpsrr);
Fps(instrumentInds,~instrumentLogicals) = DeltarInv*DeltaFpsrxtilde;
Fps(instrumentInds,instrumentInds) = DeltarInv*DeltaFpsrr;
% Define F^{pol} matrices
Fpol(~instrumentLogicals,~instrumentLogicals) = ...
    -ThetaInv*HtildeCr*DeltarInv*DeltaFpolrxtilde;
Fpol(~instrumentLogicals,instrumentInds) = ...
    -ThetaInv*HtildeCr*DeltarInv*DeltaFpolrr;
Fpol(instrumentInds,~instrumentLogicals) = DeltarInv*DeltaFpolrxtilde;
Fpol(instrumentInds,instrumentInds) = DeltarInv*DeltaFpolrr;

%% CREATE THE FUNCTION HANDLE TO COMPUTE ANTICIPATED SHOCK IMPACT MATRICES
% This will allow users to compute F matrices "on the fly" up to a horizon
% of choice with control over the format of the output.
FODhandle = @(H,returnIntermediateHorizons,offsetHorizonForProjection) ...
    compute_anticipated_shock_impact_matrices_for_ODPP(...
    H,returnIntermediateHorizons,offsetHorizonForProjection,beta,Omega,...
    BOD,Fps,Fpol);

end

%% HELPER FUNCTION TO VALIDATE OPTIONS
function [algorithmIsDennis,initialisationSource,InitialisationData,...
    AlgorithmOpts] = unpack_and_validate_solution_options_structure(...
    Options)
% Validates the optional solution options input.
%
% INPUTS:   
%   -> Options (optional): structure
%       - algorithmIsDennis (optional): true/false
%       - Initialisation (optional): Initialisation choice for B/V
%           - source (optional): 'setToZeros', 'useInputModel', 
%             'setToUserInput'
%           - Binit: nx*nx matrix for `setToUserInput' option
%           - VxtildextildeInit: nxtilde*nxtilde matrix for 
%             `setToUserInput' option if algorithmIsDennis = false
%       - Algorithm (optional): structure
%           - tol (optional): tolerance for convergence
%           - maxIter (optional): maximum number of iterations
%
% OUTPUTS:  
%   -> algorithmIsDennis: true/false
%   -> initialisationSource: 'setToZeros','useInputModel','setToUserInput'
%   -> InitialisationData: structure containing initialisation data for
%      'setToUserInput':
%           - Binit: nx*nx matrix
%           - VxtildextildeInit: nxtilde*nxtilde matrix required if
%             algorithmIsDennis = false
%   -> AlgorithmOpts: 
%           - tol (optional): tolerance for convergence
%           - maxIter (optional): maximum number of iterations

%% ALGORITHM IS DENNIS
algorithmIsDennis = true;
if isfield(Options,'algorithmIsDennis')
    algorithmIsDennis = Options.algorithmIsDennis;
    if ~is_logical_scalar(algorithmIsDennis)
        error(['algorithmIsDennis field of Options structure must ',...
            'be a true/false logical scalar']);
    end
end

%% INITIALISATION SOURCE
if algorithmIsDennis
    initialisationSource = 'useInputModel';
else
    initialisationSource = 'setToZeros';
end
InitialisationData = struct;
validInitialisationSources = {...
    'setToZeros' 'useInputModel' 'setToUserInput'};
if isfield(Options,'Initialisation')
    InitialisationOpts = Options.Initialisation;
    if ~isstruct(InitialisationOpts)
        error(['Initialisation field of Options structure input must ',...
            'be a structure']);
    end
    if isfield(InitialisationOpts,'source')       
        initialisationSource = InitialisationOpts.source;
        if ~ischar(initialisationSource) && ~any(...
                strcmp(initialisationSource,validInitialisationSources))
            error(['Initialisation.source field of Options structure ',...
                'must be either ''setToZeros'', ''useInputModel'', or ',...
                '''setToUserInput''']);
        end
        InitialisationData = InitialisationOpts;
        InitialisationData = rmfield(InitialisationData,'source');
        if strcmp(initialisationSource,'setToUserInput')
            if ~isfield(InitialisationData,'Binit')
                error(['Initialisation.source field of Options ',...
                    'structure input is set to ''setToUserInput'' in ',...
                    'which case the Initialisation field of the ',...
                    'Options structure input must also contain a ',...
                    'Binit field containing an initialisation matrix'])
            end
            if ~algorithmIsDennis
                if ~isfield(InitialisationData,'Vinit')
                    error(['Initialisation.source field of Options ',...
                        'structure input is set to ''setToUserInput'' ',...
                        'and algorithm is set to VFI ',...
                        '(algorithmIsDennis=false)in which case the ',...
                        'Initialisation field of the Options ',...
                        'structure input must also contain a Vinit ',...
                        'field containing an initialisation matrix'])
                end
            end
        end
    end
end

%% ALGORITHM OPTIONS
AlgorithmOpts = struct;
if isfield(Options,'Algorithm')
    AlgorithmOpts = Options.Algorithm;
    if ~isstruct(AlgorithmOpts)
        error(['Algorithm field of Options structure input must ',...
            'be a structure']);
    end
end

end

%% HELPER FUNCTION TO INITIALISE
function [B0,V0xtildextilde] = ...
    validate_and_unpack_user_input_initialisation(...
    InitialisationData,algorithmIsDennis,nx,nxtilde)
% Creates initialisation/guess for solution.
%
% INPUTS:   
%   -> InitialisationData
%      - Binit: nx*nx matrix
%      - Vinit: nxtilde*nxtilde matrix
%   -> algorithmIsDennis: true/false
%   -> nx: number of endogenous variables
%   -> nxtilde: number of non-instrument endogenous variables
%
% OUTPUTS:  
%   -> B0: nx*nx matrix initialisation for law of motion
%   -> V0xtildextilde: nxtilde*nxtilde initialisation for value function  

%% UNPACK INITIALISATION FOR LAW OF MOTION
B0 = InitialisationData.Binit;
if ~is_finite_real_two_dimensional_numeric_matrix(B0) || ...
        size(B0,1)~=nx || size(B0,2)~=nx
    error(['Initialisation.Binit field of Options ',...
        'structure input must contain a finite real ',...
        'two-dimensional matrix of the same dimension ',...
        'as that in the model'])
end

%% UNPACK INITIALISATION FOR VALUE FUNCTION
V0xtildextilde = [];
if ~algorithmIsDennis
    V0xtildextilde = InitialisationData.Vinit;
    if ~is_finite_real_two_dimensional_numeric_matrix(V0xtildextilde) ||...
            size(V0xtildextilde,1)~=nxtilde || ...
            size(V0xtildextilde,2)~=nxtilde
        error(['Initialisation.Vinit field of Options ',...
            'structure input must contain a finite real ',...
            'two-dimensional matrix of the same dimension ',...
            'as the number of endogenous variables in the ',...
            'model minus the number of instruments being ',...
            'used'])
    end
end

end