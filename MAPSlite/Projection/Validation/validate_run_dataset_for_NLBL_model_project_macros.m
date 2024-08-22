function validate_run_dataset_for_NLBL_model_project_macros(Model,RunData)
% This helper validates forecast run datasets for NLBL model project 
% macros. It validates that the forecast run dataset structure input 
% contains all of the information required by MAPS' family of non-linear
% backward looking model projection macros given the model being used.
%
% INPUTS:
%   -> Model: valid MAPS model
%   -> RunData: forecast run dataset sructure
%
% OUTPUTS:
%   -> none
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> unpack_model
%   -> get_NLBL_model_forecast_run_dataset_configs
%   -> validate_MAPS_dataset_content
%
% DETAILS:
%   -> This macro helper is part of a set of MAPS macro helpers that
%      validate either the model and/or the dataset being used. The aim of 
%      these helpers is to avoid repetition of checking in the macros and 
%      to remove unnecessary "noise" in the code that may be a barrier to 
%      understanding what they do.
%   -> It first validates that the model input is indeed a non-linear 
%      backward looking.
%   -> It then defines the data series that must form part of an NLBL model
%      forecast dataset in order to run any of the NLBL model project 
%      macros (including imposing judgement).
%   -> It then calls a generic MAPS dataset validation routine which will
%      validate any dataset against a (correctly formatted) dataset 
%      configuration cell array given the model being used. This function
%      will throw an exception if any data that should be present in the 
%      run dataset for operation of MAPS project macros is missing.
%
% NOTES:
%   -> See <> for more details of MAPS forecast macro helpers and data 
%      validation.
%
% This version: 06/06/2011
% Author(s): Alex Haberis and Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of inputs is as expected. All inputs are
% compulsory.
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);   
elseif ~isstruct(RunData)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
end

%% VALIDATE INPUT MODEL IS NON-LINEAR
% Validate that the input model is non-linear. Throw an exception 
% if it is not.
modelIsLinearStateSpace = unpack_model(Model,{'modelIsLinearStateSpace'});
if modelIsLinearStateSpace
    errId = ['MAPS:',mfilename,':BadModelClass'];
    generate_and_throw_MAPS_exception(errId);
end

%% GET LSS MODEL FORECAST RUN DATASET CONFIGURATION
% Get the LSS model forecast run dataset configuration information from the
% relevant config file.
runDataConfig = get_NLBL_model_forecast_run_dataset_configs;

%% FIND PAST DATA INFO IN THE CONFIG
% None of the MAPS NLBL model projection macros require past data. Compute 
% the index logicals of the past data in the configuraion.
pastConfigLogicals = strcmp('Past',runDataConfig(:,1));

%% FIND ENDOGENOUS VARIABLE FORECAST INFO IN THE CONFIG
% None of the MAPS NLBL model projection macros require endogenous 
% variables over the forecats. Compute the index logicals of the past data
% in the configuraion.
forecastEndogVarConfigLogicals = (strcmp('Forecast',runDataConfig(:,1))&...
    (strcmp('NonIdentity',runDataConfig(:,2))|...
    strcmp('Identity',runDataConfig(:,2))));

%% REMOVE THE INFO FROM THE CONFIG
% Remove the past data and endogenous variable forecast data information
% from the configuration.
configInfoToRemove = (pastConfigLogicals|forecastEndogVarConfigLogicals);
runDataToValidateConfig = runDataConfig(~configInfoToRemove,:);

%% VALIDATE THE FORECAST RUN DATASET
% Call a MAPS EASE interface helper to validate that the forecast data
% contains the required information given the configuration information
% that was removed above. this helper will throw an exception if compulsory
% information is missing from the model. If it does, add the exception as
% cause to another exception which will give context for the validation
% failure.
try
    validate_MAPS_dataset_content(Model,RunData,runDataToValidateConfig);
catch DataValidationE
    errId = ['MAPS:',mfilename,':BadRunData'];
    generate_MAPS_exception_add_cause_and_throw(DataValidationE,errId);
end

end