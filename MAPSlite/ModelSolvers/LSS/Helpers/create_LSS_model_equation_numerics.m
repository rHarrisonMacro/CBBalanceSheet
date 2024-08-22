function [HB,HC,HF,PSI] = create_LSS_model_equation_numerics(...
    thetaAndSS,HBfunHandle,HCfunHandle,HFfunHandle,PSIfunHandle)            
% This function computes an LSS model's structural numeric matrices.
% It evaluates the symbolic structural matrices associated with a linear
% state space (LSS) model's model equations given a set of parameters
% and steady state expressions.
%
% INPUTS:   
%   -> thetaAndSS: cell array containing a vector of numeric parameter
%      values and steady state values (model dependent)
%   -> HBfunHandle: function handle symbolics of loadings on lagged model
%      variables
%   -> HCfunHandle: function handle symbolics of loadings on
%      contemporaneous model variables
%   -> HFfunHandle: function handle symbolics of loadings on future
%      expected model variables
%   -> PSIfunHandle: function handle symbolics of loadings on shocks
%
% OUTPUTS:  
%   -> HB: numeric matrix of loadings on lagged model variables
%   -> HC: numeric matrix of loadings on contemporaneous model variables
%   -> HF: numeric matrix of loadings on forward model variables
%   -> PSI: numeric matrix of loadings on shocks
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> evaluate_symbolic_matrices
%
% DETAILS:  
%   -> This function converts a set of LSS model equation structural
%      symbolic matrices to numeric matrices given numeric values for the
%      parameters and steady states (if applicable) of the model.
%   -> The input and output matrices relate to the following equation:
%      HB*x{t-1}+HC*x{t}+HF*x{t+1}=PSI*z{t}
%   -> The input symbolic versions of these matrices are represented in
%      MAPS using function handles.
%
% NOTES:
%   -> See <> for a description of symbolic MAPS.
%   -> Note that all input checking is left to the utility function called
%      below (with the aim of maximising peformance given the
%      implementation). 
%   -> Note also that this function does not check that the function
%      handles input and the numeric matrices output are consistent with 
%      each other or the model being used.
%
% This version: 16/05/2011
% Author(s): Alex Haberis and Matt Waldron

%% CHECK INPUTS
% Check that the number of inputs is as expected. All other input checking
% is left to the function call below.
if nargin < 5
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~is_row_cell_array(thetaAndSS) || ...
        (size(thetaAndSS,2)~=1&&size(thetaAndSS,2)~=2)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
end

%% SETUP NAMES FOR FUNCTION HANDLE INPUTS
% Setup function handle argument names for the parameters & steady state
% values depending on whether steady states are part of the input or not.
if size(thetaAndSS,2) > 1
    thetaAndSSnames = {'theta','ss'};
else
    thetaAndSSnames = {'theta'};
end

%% EVALUATE FUNCTION HANDLES
% Call the symbolic MAPS generic symbolic matrix evaluation function to
% compute the numeric matrices.
[HB,HC,HF,PSI] = evaluate_symbolic_matrices(thetaAndSS,thetaAndSSnames,...
    HBfunHandle,HCfunHandle,HFfunHandle,PSIfunHandle);

end