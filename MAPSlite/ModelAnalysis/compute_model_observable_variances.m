function [YmodelVar,YdataVar] = compute_model_observable_variances(Model,Y)
% This model analysis function computes variances of model observables.
% It computes the variances of the observables in the model and the 
% equivalent variances in the data.
% 
% INPUTS:
%   -> Model: LSS model structure
%   -> Y: nY*T matrix of model observables over a particular sample
%
% OUTPUTS:
%   -> YmodelVar: nY*1 vector of observable variances from model
%   -> YdataVar: nY*1 vector of observable variances from data
%
% DETAILS:
%   -> This model analysis function computes the variances of a model's
%      observables as implied by the structure of the model and the 
%      equivalent set of variances implied by data for those observables.
%   -> It uses a more general function to compute auto-covariances and then
%      filters the variances out those outputs. This general function also
%      contains error handling for the inputs.
%
% NOTES:
%   -> Please see the MAPS user guide for more information on model
%      analysis functionality in MAPS.
%
% This version: 09/01/2014
% Author(s): Matt Waldron

%% CHECK INPUTS
% Note that checking of the input data types is left to the more general
% function called below.
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId);
end

%% CALL UNDERLYING MODEL ANALYSIS FUNCTION
[~,YYmodel,YYdata] = compute_model_observable_autocovariances(Model,Y,0);

%% TAKE DIAGONALS OF THE COVARIANCE MATRIX
YmodelVar = diag(YYmodel);
YdataVar = diag(YYdata);

end