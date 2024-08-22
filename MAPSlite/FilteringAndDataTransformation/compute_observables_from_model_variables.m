function Y = compute_observables_from_model_variables(x,D,G,V,w)
% This module computes data for model observables given model variables. 
% It uses the model matrices associated with the measurement equation in a 
% linear state space (LSS) model to compute a dataset of model observables 
% given a dataset of model variables (and, optionally (dependent on the 
% model), a set of measurement errors).
%
% INPUTS:   
%   -> x: nx*H matrix of model variables 
%   -> D: nY*1 vector of measurement equation constants
%   -> G: nY*nx matrix of loadings on the model variables
%   -> V (model dependent/optional): nY*nw matrix of loadings on 
%      measurement errors
%   -> w (model dependent/optional): nw*H matrix of measurement errors 
%
% OUTPUTS:  
%   -> Y: matrix of model observable data
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> is_finite_real_two_dimensional_numeric_matrix
%   -> validate_measurement_equation_matrices
%
% DETAILS:  
%   -> This forecast module computes a dataset of model observables given
%      a dataset of model variables and matrices from the linear state 
%      space model's measurement equations. 
%   -> Optionally (and dependent on the model), loadings on measurement
%      errors and a dataset of measurement errors may also be provided. 
%   -> In an arbitrary period t, the model observables are computed as 
%      Y(:,t) = D+G*x(:,t)+V*w(:,t) with measurement error and
%      Y(:,t) = D+G*x(:,t) without measurement error.
%
% NOTES:   
%   -> The purpose of this module is to avoid repetition of the same code 
%      in multiple places (eg. computing model observables in a projection 
%      and in impulse response functions).
%   -> See <> for a description of MAPS forecast modules and macros.
%
% This version: 22/02/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of the inputs is as expected. This 
% function must be passed either 3 or 5 inputs.
if nargin<5 && nargin~=3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});     
elseif ~is_finite_real_two_dimensional_numeric_matrix(x)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);  
elseif nargin>4 && ~is_finite_real_two_dimensional_numeric_matrix(w)
    errId = ['MAPS:',mfilename,':BadInput5'];
    generate_and_throw_MAPS_exception(errId);
end

%% VALIDATE MEASUREMENT EQUATION MATRICES
% Call a MAPS module helper to validate that the measurement equation 
% matrices are as expected and consistent with one another. The input
% arguments passed depend on whether or not the inputs to this function
% include measurement errors.
if nargin > 3
    validate_measurement_equation_matrices(D,G,V);
else
    validate_measurement_equation_matrices(D,G);
end

%% CHECK CONSISTENCY OF INPUTS WITH EACH OTHER
% Check that the dimensions of the matrices passed in are consistent with
% each other and the data.
T = size(x,2);
if size(x,1) ~= size(G,2)
    errId = ['MAPS:',mfilename,':IncompatibleGx'];
    generate_and_throw_MAPS_exception(errId);        
elseif nargin > 4
    if size(w,1) ~= size(V,2)
        errId = ['MAPS:',mfilename,':IncompatibleVw'];
        generate_and_throw_MAPS_exception(errId);
    elseif size(w,2) ~= T
        errId = ['MAPS:',mfilename,':Incompatiblexw'];
        generate_and_throw_MAPS_exception(errId);
    end
end

%% COMPUTE MODEL OBSERVABLES
% Compute the model observables either as a function of model variables or,
% if 5 inputs were passed in, a function of model variables and measurement
% errors.
Y = D*ones(1,T)+G*x;
if nargin == 5
    Y = Y+V*w;
end

end