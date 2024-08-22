function ValidationE = validate_incomplete_model_dataset(...
    modelDataset,modelMetadata,ValidationE)
% This helper validates that a dataset is real & does not contain infs.
% If any non-real or inf (and -inf) are found, then it throws an exception 
% which details the cause of the validation failure and the offending 
% series' metadata.
%
% INPUTS:
%   -> modelDataset: matrix of data for a particular model variable type   
%   -> modelMetadata: cell string array of model metadata
%   -> ValidationE (optional): existing exception to add causes to
%
% OUTPUTS:  
%   -> ValidationE: new or updated exception
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> validate_model_dataset_is_real
%   -> validate_model_dataset_does_not_contain_inf
%
% DETAILS:  
%   -> This forecast macro helper validates that an incomplete dataset for 
%      a particular variable type in the model (eg model variables) does 
%      not contain non-real or inf (and -inf) values.
%   -> The dataset is assumed to be potentially incomplete, meaning that
%      missing values may be represented by NaNs. In addition, the dataset
%      need not include all the variables of a given type (just those 
%      described in the meatdata passed in).
%   -> This ia a useful check to make after the transformation of an 
%      incomplete set of data from one space to another. For example,
%      transforming observables from raw observable space to model
%      observable space during the course of imposing judgement on a 
%      forecast.
%   -> If all the data passed in is real and non-inf (and the number/shape 
%      of the inputs is as expected), then this function returns no output.
%   -> The validation is done in two stages. First, this function calls a 
%      routine which validates that the dataset is real. Second, it calls a 
%      routine which validates that the dataset does not contain inf (or 
%      -inf) values.    
%   -> If some of the data fails an exception, then this function adds the
%      causes of failure to an exception. That exception is either passed
%      in as input (4 inputs must be specified) or constructed (3 inputs
%      must be specified).
%   -> Note that this function returns an exception as output regardless of
%      whether the data passed or failed validation.
%
% NOTES:   
%   -> See <> for a description of data errors exception handling in MAPS.
%
% This version: 22/02/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of the inputs is as expected. The 3rd
% input is optional.
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~iscellstr(modelMetadata) || size(modelMetadata,2)~=1
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif nargin==3 && ~strcmp(class(ValidationE),'MException')
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
end

%% CONSTRUCT AN EXCEPTION TO ADD CAUSES TO IF NONE WERE PASSED IN
% If only three inputs were passed in, construct an exception to add causes
% to as encountered below.
if nargin < 3
    errId = ['MAPS:',mfilename,':ValidationFailure'];
    ValidationE = generate_MAPS_exception(errId);
end

%% VALIDATE BASIC SHAPE OF THE DATASET
% Validate that the dataset is a two-dimensional matrix.
if ~isnumeric(modelDataset) || ndims(modelDataset)~=2
    errId = ['MAPS:',mfilename,':BadDatasetFormat'];
    generate_MAPS_exception_add_as_cause_and_throw(ValidationE,errId);
end

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

end