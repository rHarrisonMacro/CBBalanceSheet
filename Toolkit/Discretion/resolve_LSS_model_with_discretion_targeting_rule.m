function [HBOD,HCOD,HFOD,PSIOD,BOD,PHIOD,FOD,ss] = ...
    resolve_LSS_model_with_discretion_targeting_rule(ModelwTgtRule,...
    paramVec,ODinfo,Options)
% Resolves MAPS LSS model with a discretion targeting rule.

% INPUTS:   
%   -> ModelwTgtRule: MAPS LSS model with a discretion targeting rule
%   -> paremVec: parameter vector ro resolve the model with
%   -> ODinfo: structure of info required to solve the OD problem
%       - policyEqNames: single string or column cell string array of
%         targeting rule names
%       - policyShockMnems: '' (reflecting that policy shocks have been
%         removed)
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
%   -> HBOD: loadings on lagged variables in structural eqs
%   -> HCOD: loadings on contemporaneous variables in structural eqs
%   -> HFOD: loadings on one-period ahead variables in structural eqs
%   -> PSIOD: loadings on shocks in structural eqs
%   -> BOD: state transition for model solved under optimal discretion
%   -> PHIOD: loadings on shocks in optimal discretion solution
%   -> FOD: forward loading matrix set to zeros to reflect that the
%      standard solution for the forward loadings matrix is not valid
%   -> ss: updated vector of numerically "solved out" steady state values
%      (set to [] if the model does not have steady state equations) 
%
% DETAILS: 
%   -> This function is designed for use as a custom solver.  It permits a
%      model with a discretion targeting rule to be used with the rest of 
%      MAPS (with some extensions). 
%   -> The discretion targeting rule is computed using a MAPS 
%      implementation of the Dennis (2007) algorithm.
%   -> The format of the model and associated limitations is the same as in
%      "create_LSS_model_with_discretion_targeting_rule".
%
% NOTES: 
%   -> Error handling is kept to a bare minimum, reflecting the designed
%      usage as a custom solver, whereby the inputs are "pre-checked".
%
% This version: 05/03/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 3
    error('This function requires 3 or 4 inputs')
end

%% HANDLE OPTIONAL ALGORITHM OPTIONS INPUT
if nargin < 4
    Options = struct;
end

%% SOLVE OUT STEADY STATE IF NECESSARY
modelHasSteadyStateEqs = unpack_model(...
    ModelwTgtRule,{'modelHasSteadyStateEqs'});
if modelHasSteadyStateEqs
    SSfunHandle = unpack_model(ModelwTgtRule,{'SSfunHandle'});
    ss = solve_LSS_model_steady_state(SSfunHandle,paramVec);
    paramsAndSS = {paramVec,ss};
else
    ss = [];
    paramsAndSS = {paramVec};
end

%% UNPACK, EVALUATE & REPACK FUNCTION HANDLES FOR STRUCTURAL EQUATIONS
[HBfunHandle,HCfunHandle,HFfunHandle,PSIfunHandle] = unpack_model(...
    ModelwTgtRule,...
    {'HBfunHandle','HCfunHandle','HFfunHandle','PSIfunHandle'});
[HB,HC,HF,PSI] = create_LSS_model_equation_numerics(paramsAndSS,...
    HBfunHandle,HCfunHandle,HFfunHandle,PSIfunHandle);
ModelwTgtRule = pack_model(...
    ModelwTgtRule,{'HB','HC','HF','PSI'},{HB,HC,HF,PSI});

%% SOLVE MODEL UNDER DISCRETION
[BOD,PHIOD,HBOD,HCOD,HFOD,PSIOD] = ...
    solve_LSS_model_under_optimal_discretion(ModelwTgtRule,ODinfo,Options);

%% SET FORWARD LOADINGS MATRIX AS ZEROS
nx = size(BOD,1);
FOD = zeros(nx,nx);

end