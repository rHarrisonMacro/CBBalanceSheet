function nx = interrogate_LSS_model_variable_initial_conditions(x0,B)
% This helper interrogates LSS model variable initial conditions. 
% It examines the linear state space (LSS) model variable initial 
% conditions input to valdiate whether or not they are valid for use in
% MAPS LSS model forecast modules. It also returns the dimensions of the
% model variables for use in MAPS modules as required.
%
% INPUTS:
%   -> x0: vector of model variable initial conditions
%   -> B: matrix of loadings on lagged model variables from the solution
%
% OUTPUTS
%   -> nx: dimensions of the model variables
%
% CALLS: 
%   -> generate_and_throw_MAPS_exception
%   -> is_data_finite_real_column vector
%
% DETAILS: 
%   -> This helper validates model variable initial conditions for use in 
%      MAPS forecast modules. It checks that the data is valid for use in
%      MAPS and that the dimensions of the initial conditions are
%      consistent with the dimensions of the model variables in the model.
%
% NOTES:   
%   -> See <> for a description of MAPS forecast helpers and data 
%      validation.
%
% This version: 14/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number of inputs is as expected. Check that the model
% variable initial conditions are a column vector of real, finite (non-NaN
% and non-inf) data.
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});  
elseif ~is_finite_real_numeric_column_vector(x0)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
end

%% COMPUTE DIMENSIONS OF THE MODEL VARIABLES IN THE MODEL
% Compute the dimensions of the model variables in the model as the number 
% of rows in the lagged model variables matrix from the solution, where:
% x{t} = B*x{t-1}+PHI*z{t}  
nx = size(B,1);

%% CHECK DIMENSIONS OF THE MODEL VARIABLES ARE CONSISTENT
% Check that the dimensions of the model variable initial conditions input
% are consistent with the dimensions of the model variables in the model.
if size(x0,1) ~= nx
    errId = ['MAPS:',mfilename,':IncompatibleBx0'];
    generate_and_throw_MAPS_exception(errId);
end

end