function etatfNew = ...
    update_time_trends_forecast_dataset_with_judgements(...
    Model,etatf,etatfjCell)
% This helper updates a LSS model trends forecast dataset for judgements.
% It takes the values in the judgements input and overlays them on to the
% existing linear state space (LSS) model deterministic time-varying trends 
% forecast dataset.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> etatf: time-varying trends forecast dataset
%   -> etatfjCell: cell array of judgments to time-varying trends
%
% OUTPUTS:  
%   -> etatfNew: updated time-varying trends forecast dataset
%
% CALLS:    
%   -> generate_and_throw_MAPS_exception
%   -> unpack_model
%   -> update_forecast_dataset_with_judgements
%
% DETAILS:  
%   -> This helper updates a linear state space model deterministic time-
%      varying trends forecast dataset for judgements.
%   -> It updates the forecast data values for each time-varying trend over 
%      which judgement is being applied.
%           
% NOTES:    
%   -> Please see <> for a description of MAPS modules.
%   -> See also <> for a description of imposing judgement in MAPS.
%   -> This function leaves most error checking to the underlying function
%      called (to avoid repetition).
%           
% This version: 04/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of the inputs is as expected. All inputs
% are compulsory.
if nargin < 3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);    
end
    
%% UNPACK SHOCKS MEATDATA FOR JUDGEMENTS
% Unpack the model meatdata for the shocks that is used in the judgement
% structure to specify which shock(s) the judgement(s) is to be applied to.
etatMnems = unpack_model(Model,{'etatMnems'});

%% UPDATE SHOCKS
% Update each shock type over which judgement is being made. Catch any
% exceptions and add an appropriate message detailing which type of
% judgement could not be applied.
try
    etatfNew = update_forecast_dataset_with_judgements(etatMnems,...
        etatf,etatfjCell);
catch DatasetUpdateE
    errId = ['MAPS:',mfilename,':TimeTrendsDatasetUpdateFailure'];
    generate_MAPS_exception_add_cause_and_throw(...
        DatasetUpdateE,errId);
end
    

end