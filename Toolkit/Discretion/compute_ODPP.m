function [xf,ODmodel,diagnostics] = ...
    compute_ODPP(Model,ODinfo,x0,Shocks,Options)
% Computes an Optimal Discretionary Policy Projection.
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
%       - Constraints: (optional/application specific) sructure describing
%         S*r >= b:
%           - instrumentMnems: single string or column cell string array
%           - instrumentCoeffs: numeric column vector of coeffs on 
%           - constants: numeric column vector of constraint constants
%   -> x0: nx*1 vector of initial conditions for model variables
%   -> Shocks: struture with at least one of 
%       - anticipated (model dependent): nz*H matrix of shocks
%       - unanticipated: nz*H matrix of shocks (can be zeros for non-
%         judgemental projection)
%   -> Options (optional): structure
%       - Solution (optional): structure
%           - algorithmIsDennis (optional): true/false
%           - Initialisation (optional): Initialisation choice for B/V
%               - source (optional): 'setToZeros', 'useInputModel', 
%               'setToUserInput'
%               - Binit: nx*nx matrix for `setToUserInput' option
%               - Vinit: nxtilde*nxtilde matrix for 
%               `setToUserInput' option if algorithmIsDennis = false
%           - Algorithm (optional): structure
%               - tol (optional): tolerance for convergence
%               - maxIter (optional): maximum number of iterations
%       - Projection (optional): structure (NB. relevant for constrained
%         projection only)
%           - initialisationFuncHandle (optional): function to use to 
%             initialise constraint binding indicator function for BPY
%           - BPY (optional): options for the BPY algorithm
%               - tol: tolerance for convergence
%               - maxIter: maximum number before giving up
%               - reportProgress: false/true
%               - useBestOnNonConvergence: false/true
%
% OUTPUTS
%   -> xf: nx*H matrix of projected model variables
%   -> ODmodel: information about the OD problem and resulting solution --
%      ODinfo structure with the following fields added:
%       - policyEqLogicals: vector of logicals that locate policy eqs
%       - policyEqInds: vector of indices that locate policy eqs
%       - policyShockLogicals: vector of logicals locating policy shocks
%       - policyShockInds: vector of indices locating policy shocks
%       - instrumentLogicals: vector of logicals locating the instruments
%       - instrumentInds: vector of indices that locate the instruments
%       - Q: instruments weights matrix
%       - W: target weights matrix
%       - B: state transition for model solved under optimal discretion
%       - PHI: loadings on shocks in optimal discretion solution
%       - FODhandle: function handle that can be used to compute the
%         horizon-varying anticipated shock loadings for a user choice of
%         horizon (and treatment of the intervening periods and starting 
%         point -- i.e. does counting begin at F_1 or F_0) 
%       - Vxtildextilde: quadratic valud function coefficients on non-
%         instrument endogenous variables
%
% DETAILS:
%   -> Computes an Optimal discretionary policy projection (ODPP) 
%      in the LQ framework.
%   -> The overarching logic is the same as in Svennson and Tetlow -- the
%      optimal discretionary policy projection is conditional on the 
%      initial conditions, x0, and the judgements implicit in the shocks.
%   -> The algorithm proceeds in two steps.  First, the model is solved for
%      optimal discretionafry policy using either the Dennis algorithm
%      (which is the default) or VFI (which can be selected via
%      Options.Solution).  Second, a projection is computed using one of
%      the two algorithms, depending on whether or not there are constraints 
%      on the instruments (in ODinfo.Constraints).  Note that if there are 
%      instrument constraints, then there are some optional projection 
%      options that can be passed in as Options.Projection.  Note also that 
%      there is an assumption that the constraints do not bind in steady 
%      state.
%
% NOTES:   
%   -> Note that error handling is kept to a minimum (with some (but only) 
%      some of it delegated to the functions called below).
%
% This version: 05/03/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 4
    error('This function requires 4 or 5 inputs')
end

%% HANDLE OPTIONAL OPTIONS INPUT
SolOpts = struct;
ProjectOpts = struct;
if nargin > 4
    check_field_names_in_structure(...
        Options,'optional options input','',{'Solution';'Projection'});
    if isfield(Options,'Solution')
        SolOpts = Options.Solution;
    end
    if isfield(Options,'Projection')
        ProjectOpts = Options.Projection;
    end        
end

%% SOLVE MODEL WITH OPTIMAL DISCRETIONARY POLICY
[BOD,PHIOD,HBOD,HCOD,HFOD,PSIOD,FODhandle,ODinfo,Vxtildextilde] = ...
    solve_LSS_model_under_optimal_discretion(Model,ODinfo,SolOpts);

%% DETERMINE TYPE OF PROJECTION
if isfield(ODinfo,'Constraints')
    projectionIsSubjectToConstraints = true;
else
    projectionIsSubjectToConstraints = false;
end

%% COMPUTE AN UNCONSTRAINED PROJECTION
if ~projectionIsSubjectToConstraints
    xf = project_using_OD_model(BOD,PHIOD,FODhandle,x0,Shocks);   
end

%% COMPUTE A CONSTRAINED PROJECTION
% Note that Htail is the number of periods to use to compute expectations
% at the end of the projection horizon.  There is a presumption that there
% will have been convergence to a regime in which there are no binding
% instrument constraints within the H period window.
if projectionIsSubjectToConstraints
    Htail = 0;
    [xf,Jf] = project_using_OD_model_subject_to_instrument_bound_constraints(...
        ODinfo,HBOD,HCOD,HFOD,PSIOD,BOD,Vxtildextilde,x0,Shocks,Htail,...
        ProjectOpts);
end

%% CONSTRUCT OD MODEL OUTPUT
ODmodel = ODinfo;
ODmodel.B = BOD;
ODmodel.PHI = PHIOD;
ODmodel.FODhandle = FODhandle;
ODmodel.Vxtildextilde = Vxtildextilde;

%% ADD DIAGNOSTICS
diagnostics = struct;
if projectionIsSubjectToConstraints
    diagnostics.Jf = Jf;
end

end