function validate_MAPS_dataset_content(Model,MAPSdata,MAPSdataConfig)
% This helper validates the content of a generic MAPS dataset structure.
% It is used as a helper in the validation of the content of MAPS datasets 
% in EASE interfaces and MAPS macros.
%
% INPUTS:
%   -> Model: MAPS model structure
%   -> MAPSdata: MAPS dataset structure
%   -> MAPSdataConfig: configuration information for the dataset structure
%
% OUTPUTS:
%   -> none
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> unpack_model
%   -> unpack_data_from_MAPS_dataset
%   -> generate_MAPS_exception_and_add_as_cause
%
% DETAILS:
%   -> This helper validates the content of a generic MAPS dataset as
%      described by the input configuration information.
%   -> It searches the dataset input for all possible data series (as
%      defined in the configuration information) and attempts to unpack
%      them one by one if the configuration information states that they 
%      are compulsory components of the dataset or if they are compulsory 
%      given the model being used.
%   -> If the unpack data helper returns empty for these data series, then
%      this macro throws an exception detailing which of the components of
%      the dataset input were missing.
%
% NOTES:
%   -> See xxxxxx for a discussion of MAPS datasets and their validation.
%
% This version: 16/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number & shape of inputs is as expected. All inpurs are 
% compulsory. The EASE data and configuration must be cell arrays and the 
% MAPS data configuration must be a four-column cell array. 
if nargin < 3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);      
elseif ~isstruct(MAPSdata)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif ~iscellstr(MAPSdataConfig) || ndims(MAPSdataConfig)~=2 || ...
        size(MAPSdataConfig,2)<7
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
end

%% COMPUTE THE NUMBER OP POSSIBLE DATA TYPES IN THE MAPS DATASET
% Compute the number of possible data types in the MAPS dataset as the
% number of separate rows in the MAPS data configuration (each of which
% describes the configuration for a particular type).
nPossDataTypes = size(MAPSdataConfig,1);

%% SETUP A MASTER EXCEPTION
% Setup a master exception to add causes to as encountered below.
masterErrId = ['MAPS:',mfilename,':DataMissing'];
DataMissingE = generate_MAPS_exception(masterErrId);

%% CHECK THE CONTENT OF THE MAPS DATASET FOR EACH POSSIBLE DATA TYPE
% Run through each possible data type, checking whether or not that type is
% either outright compulsory or compulsory because the model 
% characteristics condition is met. If the data is compulsory, attempt to
% unpack it. If it returns empty, then add an exception with details of
% which datset was missing and the treatment of that dataset (eg compulsory
% or compulsory if ...).
for iType = 1:nPossDataTypes
    iTypeSeriesIds = MAPSdataConfig(iType,1:3);
    iTypeContentType = MAPSdataConfig{iType,7};
    iTypeShouldExist = false;
    if strcmp(iTypeContentType,'compulsory')
        iTypeShouldExist = true;
    elseif strfind(iTypeContentType,'compulsory if')
        iTypeContentCondition = strrep(...
            iTypeContentType,'compulsory if ','');
        iTypeShouldExist = unpack_model(Model,{iTypeContentCondition});
    end
    if iTypeShouldExist
        try
            modelDataset = unpack_data_from_MAPS_dataset(...
                MAPSdata,iTypeSeriesIds);
            if isempty(modelDataset)
                errId = [masterErrId,':Instance'];
                DataMissingE = generate_MAPS_exception_and_add_as_cause(...
                    DataMissingE,errId,...
                    [iTypeSeriesIds {iTypeContentType}]);
            end
        catch DataUnpackE
            errId = ['MAPS:',mfilename,':DataUnpackFailure'];
            generate_MAPS_exception_and_add_cause_and_throw(...
                DataUnpackE,errId,iTypeSeriesIds);
        end
    end
end

%% THROW ANY EXCEPTION(S) ENCOUNTERED
% Throw the master exception created above if any data content exceptions
% were encountered.
if ~isempty(DataMissingE.cause)
    throw(DataMissingE);
end

end