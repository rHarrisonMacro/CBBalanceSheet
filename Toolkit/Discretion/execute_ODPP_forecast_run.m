function [ODPPrun,ODmodel] = execute_ODPP_forecast_run(...
    Model,BaseRun,ODinfo,Options)
% Executes an Optimal Discretionary Policy Projection forecast run.
%
% INPUTS:   
%   -> Model: MAPS LSS model
%   -> BaseRun: baseline MAPS LSS model forecast run data structure
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
%             instruments
%           - constants: numeric column vector of constraint constants
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
%   -> ODPPrun: forecast run with policy acting under discretion
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
%       - HB: loadings on lagged variables in structural eqs
%       - HC: loadings on contemporaneous variables in structural eqs
%       - HF: loadings on one-period ahead variables in structural eqs
%       - PSI: loadings on shocks in structural eqs
%       - FODhandle: function handle that can be used to compute the
%         horizon-varying anticipated shock loadings for a user choice of
%         horizon (and treatment of the intervening periods and starting 
%         point -- i.e. does counting begin at F_1 or F_0) 
%       - Vxtildextilde: quadratic valud function coefficients on non-
%         instrument endogenous variables
%
% DETAILS:
%   -> Takes a baseline forecast run and computes a forecast with policy
%      acting under optimal discretion conditional on the initial
%      conditions and judgements embodied in the baseline forecast run.
%   -> See compute_ODPP and the functions therein for further details.
%
% NOTES:   
%   -> Note that most error handling is delegated to the functions called 
%      below.
%
% This version: 06/03/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 3
    error('This function requires 3 inputs')
end

%% HANDLE OPTIONAL ALGORITHM OPTIONS INPUT
if nargin < 4
    Options = struct;
end

%% COMPUTE ODPP
xT = BaseRun.Constraint.modelVariables;
Shocks = BaseRun.Forecast.Shocks;
[xf,ODmodel] = compute_ODPP(Model,ODinfo,xT,Shocks,Options);

%% PACK MODEL VARIABLE PROJECTIONS
ODPPrun = BaseRun;
ODPPrun.Forecast.modelVariables = xf;

%% PROJECT MODEL OBSERVABLES
modelHasMeasurementEqs = unpack_model(Model,'modelHasMeasurementEqs');
if modelHasMeasurementEqs
    [D,G] = unpack_model(Model,{'D','G'});
    Yf = compute_observables_from_model_variables(xf,D,G);
    ODPPrun.Forecast.modelObservables = Yf;
end

%% PROJECT RAW OBSERVABLES
modelHasDataTransformationEqs = unpack_model(...
    Model,'modelHasDataTransformationEqs');
if modelHasDataTransformationEqs
    [RTfunHandle,modelHasTimeVaryingTrends] = unpack_model(...
        Model,{'RTfunHandle','modelHasTimeVaryingTrends'});
    YtildeT = BaseRun.Constraint.rawObservables;
    if modelHasTimeVaryingTrends
        etatT = BaseRun.Constraint.timeVaryingTrends;
        etatf = BaseRun.Forecast.timeVaryingTrends;
        Ytildef = transform_observables_from_model_to_raw_space(...
            RTfunHandle,Yf,YtildeT,etatf,etatT);
    else
        Ytildef = transform_observables_from_model_to_raw_space(...
            RTfunHandle,Yf,YtildeT);  
    end
    validate_retransformed_raw_observable_dataset(Model,Ytildef);
    ODPPrun.Forecast.rawObservables = Ytildef;
end

end