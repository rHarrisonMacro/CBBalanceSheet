function EASErunData = ...
    convert_forecast_run_dataset_from_MAPS_to_EASE_format(...
    Model,MAPSrunData)
% This helper converts forecast run data output from MAPS to EASE format.
% It should be used as a helper to all EASE interfaces in MAPS that take as
% output a forecast run dataset from MAPS. It converts the structure
% forecast run dataset from MAPS to an equivalent cell array dataset for
% EASE.
%
% INPUTS:
%   -> Model: valid MAPS model
%   -> MAPSrunData: structure of forecast run data from MAPS
%
% OUTPUTS:
%   -> EASErunData: cell array of forecast run data for EASE
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> unpack_model
%   -> get_LSS_model_forecast_run_dataset_configs
%   -> get_NLBL_model_forecast_run_dataset_configs
%   -> convert_dataset_from_MAPS_to_EASE_format
%
% DETAILS:
%   -> This helper translates the MAPS structure format for a forecast run
%      to a cell array for use in EASE.
%   -> The dataset input and the translation depends on whether it is
%      associated with a linear state space model or a non-linear
%      backward-looking model.
%   -> EASE-formatted datasets are stored as cell arrays with each
%      indiviudal variable's data having its own row in the cell array
%      (where the first few columns provide metadata about the data type;
%      the penultimate column provides metadata about the individual series
%      and the final column contains the time series of data itself).
%   -> MAPS-formatted datasets are stored as structures with the fields 
%      (and sub fields) of the structures equivalent to the series metadata 
%      in the EASE cell array. For example (and if the config information 
%      dictates), forecast data for raw observables would be stored as an 
%      nYtilda*H matrix (where H is the forecast horizon) in 
%      MAPSrunData.Forecast.rawObservables.
%
% NOTES:
%   -> See xxxxxx for a discussion of EASE- & MAPS-represented forecast run
%      dataset formats in EASE and MAPS.
%
% This version: 14/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number & shape of inputs is as expected. All inputs are
% compulsory. The EASE dataset must be a cell array.
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(MAPSrunData)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
end

%% UNPACK MODEL CLASS INFO
% Unpack the model is linear state space flag to be used below to call the
% correct dataset configuration.
modelIsLinearStateSpace = unpack_model(Model,{'modelIsLinearStateSpace'});

%% GET DATA CONFIGURATIONS
% Call the configuration function to get the EASE and MAPS data
% configurations according to the model type.
if modelIsLinearStateSpace
    [MAPSrunDataConfig,EASErunDataConfig] = ...
        get_LSS_model_forecast_run_dataset_configs;
else
    [MAPSrunDataConfig,EASErunDataConfig] = ...
        get_NLBL_model_forecast_run_dataset_configs;
end

%% CONVERT EASE PROVISIONAL RUN DATA TO MAPS FORMAT
% Pass the above configurations, model and data inputs to a helper function
% which constructs the MAPS formatted data.
EASErunData = convert_dataset_from_MAPS_to_EASE_format(...
    Model,MAPSrunData,MAPSrunDataConfig,EASErunDataConfig);

end