function validate_data_in_MAPS_dataset(Model,MAPSdata,MAPSdataConfig,T,H)
% This helper validates the data in a generic MAPS dataset structure.
% It is used as a helper in the validation of data in MAPS datasets in
% EASE interfaces and MAPS macros.
%
% INPUTS:
%   -> Model: MAPS model structure
%   -> MAPSdata: MAPS dataset structure
%   -> MAPSdataConfig: configuration information for the dataset structure
%   -> T: horizon of data over the past
%   -> H: horizon of data over the forecast
%
% OUTPUTS:
%   -> none
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> unpack_data_from_MAPS_dataset
%   -> validate_model_matrix_dataset (sub-function)
%   -> validate_model_cell_dataset (sub-function)
%
% DETAILS:
%   -> This helper validates the data in a generic MAPS dataset structure 
%      as described by the input configuration information.
%   -> It searches the dataset input for all possible data series (as
%      defined in the configuration information) and attempts to unpack
%      each of them. If the data exists in the dataset then its content is
%      validated (depending on whether the dataset should be a complete
%      model matrix dataset or an incomplete cell array).
%   -> In either case, an exception is generated and added as cause to a
%      master exception if the data contains invalid numbers (like inf) or
%      if it does not have the correct time dimension (as specified in the
%      input and configuration).
%   -> If the dataset fails validation, the master exception is thrown
%      detailing all causes of the failure. If the dataset passes
%      validation, this function returns nothing.
%
% NOTES:
%   -> See xxxxxx for a discussion of MAPS datasets and their validation.
%
% This version: 21/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of inputs is as expected. All inputs are
% compulsory.
if nargin < 5
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(MAPSdata)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId); 
elseif ~iscellstr(MAPSdataConfig) || ndims(MAPSdataConfig)~=2 || ...
        size(MAPSdataConfig,2)<8
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);  
elseif ~is_positive_real_integer(T)
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);   
elseif ~is_positive_real_integer(H)
    errId = ['MAPS:',mfilename,':BadInput5'];
    generate_and_throw_MAPS_exception(errId);
end      

%% SETUP A MASTER EXCEPTION
% Setup a master exception to add causes to as encountered below.
masterErrId = ['MAPS:',mfilename,':BadData'];
BadDataE = generate_MAPS_exception(masterErrId);

%% UNPACK & VALIDATE DATA
% Run through each of the data series to validate. Attempt to unpack them
% one-by-one. If they exist in the dataset input, call a sub-function below
% to validate them based on their type. Capture any exceptions encountered
% and add them to the master exception.
nDataToValidate = size(MAPSdataConfig,1);
for iData = 1:nDataToValidate
    iDataIds = MAPSdataConfig(iData,1:3);
    modelDataset = unpack_data_from_MAPS_dataset(MAPSdata,iDataIds);
    if ~isempty(modelDataset)
        iDataStorageType = MAPSdataConfig{iData,4};
        iDataStorageSubType = MAPSdataConfig{iData,5};
        iDataModelMetadataStr = MAPSdataConfig{iData,6};
        iDataTimeDims = eval(MAPSdataConfig{iData,8});
        errId = ['MAPS:',mfilename,':BadDataset'];
        BadDatasetE = generate_MAPS_exception(errId,iDataIds);
        try
            switch iDataStorageType
                case 'matrix'
                    validate_model_matrix_dataset(Model,modelDataset,...
                        BadDatasetE,iDataStorageSubType,...
                        iDataModelMetadataStr,iDataTimeDims);
                case 'cell'
                    validate_model_cell_dataset(Model,modelDataset,...
                        BadDatasetE,iDataStorageSubType,...
                        iDataModelMetadataStr,iDataTimeDims);
                otherwise
                    errId = ['MAPS:',mfilename,':UnknownStorageType'];
                    generate_and_throw_MAPS_exception(...
                        errId,{iDataStorageType});
            end
        catch BadDataEi
            BadDataE = addCause(BadDataE,BadDataEi);
        end
    end
end

%% THROW ANY EXCEPTION(S) ENCOUNTERED
% Throw the master exception created above if any data exceptions were 
% encountered.
if ~isempty(BadDataE.cause)
    throw(BadDataE);
end

end

%% FUNCTION TO VALIDATE MODEL MATRIX DATASET
function validate_model_matrix_dataset(Model,modelDataset,...
    BadDatasetE,dataStorageSubType,modelMetadataStr,timeDims)
% This sub-routine validates complete model datasets in matrices.
% Complete model datasets are supposed to contain time series for all the
% variables of the specified type over a homogenous time horizon with no
% inf, imaginary or missing (NaN) values.
%
% INPUTS:
%   -> Model: MAPS model structure
%   -> MAPSdata: MAPS dataset structure
%   -> BadDatasetE: exception to add causes to as encountered below
%   -> dataStorageSubType: string identifier for data storage sub-type
%   -> modelMetadataStr: string identifier for model metadata
%   -> timeDims: expected time dimension of the data
%
% OUTPUTS:
%   -> none
%
% CALLS:
%   -> validate_complete_model_dataset
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception_add_as_cause_and_throw

%% VALIDATE DATASET CONTENT
% If the data sub-type is specified as "numeric", call another MAPS helper
% to validate that the dataset contains only real, finite, non-NaN numbers.
% This function will also validate that the dimensions of the data are
% consistent with the dimensions of the model metadata associated with the
% data. If the data sub-type is not numeric, throw an exception for an
% unhandled sub-type.
switch dataStorageSubType
    case 'numeric'
        validate_complete_model_dataset(...
            Model,modelDataset,modelMetadataStr,BadDatasetE);
    otherwise
        errId = ['MAPS:',mfilename,':UnknownStorageSubType'];
        generate_and_throw_MAPS_exception(errId,{dataStorageSubType});
end

%% VALIDATE DATASET TIME DIMENSIONS
% Check that the dataset has the same time dimension (number of columns) as 
% expected.
if size(modelDataset,2) ~= timeDims
    errId = ['MAPS:',mfilename,':BadDataDims'];
    generate_MAPS_exception_add_as_cause_and_throw(...
        BadDatasetE,errId,{num2str(timeDims)}); 
end
   
end

%% FUNCTION TO VALIDATE MODEL CELL DATASET
function validate_model_cell_dataset(Model,modelDataset,...
    BadDatasetE,dataStorageSubType,modelMetadataStr,timeDims)
% This sub-routine validates incomplete model datasets in cells.
% Incomplete cell aray model datasets are supposed to contain time series 
% for a sub-set of variables of the specified type (with metadata 
% describing which variables are included) with no inf or imaginary values 
% (missing (NaN) values are allowed).
%
% INPUTS:
%   -> Model: MAPS model structure
%   -> MAPSdata: MAPS dataset structure
%   -> BadDatasetE: Exception to add causes to as encountered below
%   -> dataStorageSubType: string identifier for data storage sub-type
%   -> modelMetadataStr: string identifier for model metadata
%   -> timeDims: expected time dimension of the data
%
% OUTPUTS:
%   -> none
%
% CALLS:
%   -> validate_cell_array_incomplete_model_dataset
%   -> validate_cell_array_model_dataset_of_logicals
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception_add_as_cause_and_throw

%% VALIDATE DATASET CONTENT
% Call another MAPS helper (depending on the sub-type of data) to validate
% that the metadata in the cell array is consistent with the model being
% used and that the data in the cell array contains only real, non-inf
% numbers. Unlike in complete MAPS matrix datasets, the data may contain
% NaNs, which are used as nulls to signify periods where there is no data
% (eg in judgments).
switch dataStorageSubType
    case 'numeric'
        validate_cell_array_incomplete_model_dataset(...
            Model,modelDataset,modelMetadataStr,BadDatasetE);
    case 'logical equivalent'
        validate_cell_array_model_dataset_of_logicals(...
            Model,modelDataset,modelMetadataStr,BadDatasetE);        
    otherwise
        errId = ['MAPS:',mfilename,':UnknownStorageSubType'];
        generate_and_throw_MAPS_exception(errId,{dataStorageSubType});
end
        

%% VALIDATE DATASET TIME DIMENSIONS
% Check that the dataset has the same time dimension (number of columns) as 
% expected.
if any(cellfun('size',modelDataset(:,2),2)-timeDims)
    errId = ['MAPS:',mfilename,':BadDataDims'];
    generate_MAPS_exception_add_as_cause_and_throw(...
        BadDatasetE,errId,{num2str(timeDims)});
end

end