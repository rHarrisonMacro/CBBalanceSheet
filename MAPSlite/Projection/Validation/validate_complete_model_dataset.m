function validate_complete_model_dataset(...
    Model,modelDataset,modelMetadataStr,ValidationE)
% This helper validates that a complete model dataset is real and finite.
% If any non-real (complex numbers with imaginary components), inf (and 
% -inf) or nan values are found, then it throws an exception which details
% the cause of the validation failure and the offending series' metadata.
%
% INPUTS:   
%   -> Model: MAPS model structure 
%   -> modelDataset: matrix of data for a particular model variable type   
%   -> modelMetadataStr: string for the variable type metadata in the model
%   -> ValidationE (optional): existing exception to add causes to
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> unpack_model
%   -> validate_model_dataset_is_real
%   -> validate_model_dataset_does_not_contain_inf
%   -> validate_model_dataset_does_not_contain_nan
%
% DETAILS:  
%   -> This forecast macro helper validates that a dataset for a particular
%      variable type in the model (eg model variables) does not contain 
%      non-real, inf (and -inf) or nan values. 
%   -> This ia a useful check to make after the transformation of a 
%      complete set of data from one space to another. For example,
%      transforming observables from raw observable space to model
%      observable space.
%   -> If all the data passed in is real and finite (and the number / shape 
%      of the inputs is as expected), then this function returns no output.
%   -> The validation is done in three stages. First, this function calls a 
%      routine which validates that the dataset is real. Second, it calls a 
%      routine which validates that the dataset does not contain inf (or 
%      -inf) values. Third, this routine extends that to test for NaN 
%      values - which should not be present in a complete dataset.    
%   -> If some of the data fails validation, then this function adds the
%      causes of failure to an exception. That exception is either passed
%      in as input (4 inputs must be specified) or constructed (3 inputs
%      must be specified).
%
% NOTES:   
%   -> See <> for a description of data errors exception handling in MAPS.
%
% This version: 22/02/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of the inputs is as expected. The 4th
% input is optional.
if nargin < 3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~ischar(modelMetadataStr)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
elseif nargin==4 && ~strcmp(class(ValidationE),'MException')
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);
end

%% CONSTRUCT AN EXCEPTION TO ADD CAUSES TO IF NONE WERE PASSED IN
% If only three inputs were passed in, construct an exception to add causes
% to as encountered below.
if nargin == 3
    errId = ['MAPS:',mfilename,':ValidationFailure'];
    ValidationE = generate_MAPS_exception(errId);
end

%% VALIDATE BASIC SHAPE OF THE DATASET
% Validate that the dataset is a two-dimensional matrix.
if ~isnumeric(modelDataset) || ndims(modelDataset)~=2
    errId = ['MAPS:',mfilename,':BadDatasetFormat'];
    generate_MAPS_exception_add_as_cause_and_throw(ValidationE,errId);
end

%% UNPACK THE MODEL METADATA FOR THE ROUTINES BELOW
% Call the unpacker function to unpack the model meatdata required for
% exception building in the routines called below.
modelMetadata = unpack_model(Model,{modelMetadataStr});

%% CHECK DATASET IS REAL
% Call a routine to check that the dataset does not contain complex numbers
% with imaginary components. Catch any exceptions found and add to the 
% master exception.
try
    validate_model_dataset_is_real(modelDataset,modelMetadata);
catch NonRealE
    ValidationE = addCause(ValidationE,NonRealE);
end

%% CHECK DATASET DOES NOT CONTAIN INFS
% Call a routine to check that the dataset does not contain inf (or -inf) 
% numbers. Catch any exceptions found and add to the master exception.
try
    validate_model_dataset_does_not_contain_inf(...
        modelDataset,modelMetadata);
catch InfE
    ValidationE = addCause(ValidationE,InfE);
end

%% CHECK DATASET DOES NOT CONTAIN NANS
% Extend the incomplete dataset routines for validation that the dataset
% does not contain nans. Catch any exceptions found and add to the master
% exception.
try
    validate_model_dataset_does_not_contain_nan(...
        modelDataset,modelMetadata);
catch NanE
    ValidationE = addCause(ValidationE,NanE);
end

%% THROW ANY EXCEPTIONS CAUGHT
% If the master exception contains any causes, throw it back out of this
% function.
if ~isempty(ValidationE.cause)
    throw(ValidationE);
end

end