function validate_provisional_run_dataset_output(Model,ProvRunData,T,H)
% This helper validates an output provisional run dataset structure.
% It validates both the content of and the data in a provisional run 
% dataset from MAPS.
%
% INPUTS:
%   -> Model: MAPS model structure
%   -> ProvRunData: provisional run dataset structure
%   -> T: expected horizon of data over the past
%   -> H: expected horizon of data over the forecast
%
% OUTPUTS:
%   -> none
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> unpack_model
%   -> get_LSS_model_provisional_run_dataset_configs
%   -> get_NLBL_model_provisional_run_dataset_configs
%   -> validate_MAPS_dataset_content
%   -> validate_data_in_MAPS_dataset
%
% DETAILS:
%   -> This helper validates a MAPS provisional run dataset structure.
%   -> MAPS provisional run datasets are stored as structures with the
%      fields (and sub fields) of the structures acting as metadata for the
%      series. For example (and if the config information dictates),
%      an anticipated judgement for a single raw observable would be stored
%      as a 1*2 cell array with the first element containing model metadata
%      for the variable in question (eg 'gdpkp') and the second element
%      containing judgemental data for that variable in a 1*H vector (where 
%      H is the forecast horizon) in
%      ProvRunData.Forecast.AnticipatedFixes.rawObservables.
%   -> This function validates the dataset in two stages. It first calls a
%      helper function to check that the provisional run dataset contains 
%      data for the correct variables given the model being used. The
%      "past" component of a provisional run dataset is treated as optional 
%      (because the forecast run may have been a judgemental forecast run 
%      with no new back data introduced). However, if the dataset does 
%      include new past data, then the content becomes non-optional.
%   -> It then calls another helper to check that the data is valid
%      (regardless of the content) in that it has the expected shape, the
%      right number of time periods and contains valid data (eg no inf
%      values etc).
%   -> If the provisional run dataset fails validation for either reason, 
%      this function will throw an exception detailing the causes of the
%      validation failure.
%
% NOTES:
%   -> See xxxxxx for a discussion of MAPS forecast run datasets and their
%      validation.
%
% This version: 21/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number & shape of inputs is as expected. All inputs are
% compulsory. The EASE dataset must be a cell array.
if nargin < 4
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(ProvRunData)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif ~is_positive_real_integer(T)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
elseif ~is_positive_real_integer(H)
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);
end

%% SETUP A MASTER EXCEPTION
% Setup a master exception to add cuases to as encountered below.
errId = ['MAPS:',mfilename,':BadProvRunData'];
BadProvRunDataE = generate_MAPS_exception(errId);

%% UNPACK MODEL CLASS INFO
% Unpack the model is linear state space flag to be used below to call the
% correct dataset configuration.
modelIsLinearStateSpace = unpack_model(Model,{'modelIsLinearStateSpace'});

%% GET DATA CONFIGURATIONS
% Call the configuration function to get the MAPS forecast run dataset
% configuration according to the model type.
if modelIsLinearStateSpace
    provRunDataConfig = get_LSS_model_provisional_run_dataset_configs;
else
    provRunDataConfig = get_NLBL_model_provisional_run_dataset_configs;
end

%% AMEND CONFIG TO CHECK CONTENT
% None of the MAPS LSS model projection macros require past data to
% operate. As a result, the past data is an optional component of a
% provisional run dataset. If past data is not in the dataset, remove the
% information about past data from the configuration to avoid invalid 
% content exceptions below.
nDataTypesToCheck = size(provRunDataConfig,1);
configInfoToRemove = false(nDataTypesToCheck,1);
if ~isfield(ProvRunData,'Past')
    pastConfigLogicals = strcmp('Past',provRunDataConfig(:,1));
    configInfoToRemove = (configInfoToRemove|pastConfigLogicals);
end
provRunDataToValidateConfig = provRunDataConfig(~configInfoToRemove,:);

%% VALIDATE THE CONTENT OF THE PROVISIONAL RUN DATASET
% Call an EASE interface helper to validate that the provisional dataset
% contains all of the compulsory information given the model being used.
% This helper will throw an exception if compulsory information is missing
% from the dataset. If it does, add the exception as cause to another
% exception which will give context for the validation failure.
try
    validate_MAPS_dataset_content(...
        Model,ProvRunData,provRunDataToValidateConfig);
catch DataValidationE
    BadProvRunDataE = addCause(BadProvRunDataE,DataValidationE);
end

%% VALIDATE THE DATA IN THE PROVISIONAL RUN DATASET
% Call another EASE interface helper to validate that the data in the
% provisional run dataset is ok to output to EASE.
try
    validate_data_in_MAPS_dataset(...
        Model,ProvRunData,provRunDataToValidateConfig,T,H);
catch DataValidationE
    BadProvRunDataE = addCause(BadProvRunDataE,DataValidationE);
end

%% THROW ANY EXCEPTIONS ENCOUNTERED
% If any exceptions were encountered above the master exception will
% contain causes, so throw it back out into the calling function.
if ~isempty(BadProvRunDataE.cause)
    throw(BadProvRunDataE);
end

end