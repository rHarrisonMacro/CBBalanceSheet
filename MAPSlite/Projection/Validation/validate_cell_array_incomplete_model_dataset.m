function validate_cell_array_incomplete_model_dataset(...
    Model,modelDataset,modelMetadataStr,ValidationE)
% This helper validates an incomplete model dataset cell array.
% If the metadata in the cell array cannot be reconciled with the model 
% being used or if any non-real or inf (and -inf) are found, then it throws 
% an exception which details the cause(s) of the validation failure.
%
% INPUTS:
%   -> Model: MAPS model structure
%   -> modelDataset: two-column cell array incomplete model dataset   
%   -> modelMetadataStr: string identifier for model metadata
%   -> ValidationE (optional): existing exception to add causes to
%
% OUTPUTS:  
%   -> ValidationE: new or updated exception
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> unpack_model_metadata_and_lookup_index_numbers
%   -> generate_MAPS_exception_and_add_cause
%   -> convert_cell_array_of_vectors_to_matrix_equivalent
%   -> validate_incomplete_model_dataset
%
% DETAILS:  
%   -> This forecast macro helper validates an incomplete cell arrray 
%      dataset for a particular variable type in the model (eg model 
%      variables) does.
%   -> A MAPS incomplete cell array dataset is a two -column cell array.
%      The first column contains metadata associating each individual data
%      series with a variable in the model. The second column contains 
%      individual numeric vectors of time series data for each of the 
%      identified variables (eg {'gdpkp' [30000 31000 32000 etc]};
%   -> Each of the individual data series is assumed to be potentially 
%      incomplete, meaning that missing values may be represented by NaNs.
%   -> If the inputs are as expected, this validation helper first checks 
%      that the metadata can be reconciled with the model input. 
%   -> It then checks that the data series do not contain inf (or -inf) 
%      values or complex numbers with imaginary components.   
%   -> If the validation fails, then this function adds the cause(s) of
%      failure to an exception which it then throws. That exception is 
%      either passed in as input (4 inputs must be specified) or 
%      constructed (3 inputs must be specified).
%
% NOTES:   
%   -> See <> for a description of data errors exception handling in MAPS.
%
% This version: 21/03/2011
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
elseif nargin==3 && ~strcmp(class(ValidationE),'MException')
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

%% CHECK CELL ARRAY DATATSET CONTENT
% Check that the dataset is indeed a two-column cell array. If it is check 
% that the first column of the dataset cell array is a cell string array
% of metadata for each of the series. Check also that the second column is 
% a column cell array of numeric row vectors of time series data for each 
% of the variables identified with the metadata in the first column.
if ~iscell(modelDataset) || ndims(modelDataset)~=2 || ...
        size(modelDataset,2)~=2        
    errId = ['MAPS:',mfilename,':BadDatasetFormat'];
    generate_MAPS_exception_add_as_cause_and_throw(ValidationE,errId);
end
if ~iscellstr(modelDataset(:,1))
    errId = ['MAPS:',mfilename,':BadColumn1'];
    generate_MAPS_exception_and_add_as_cause(ValidationE,errId);
elseif ~all(cellfun(@isnumeric,modelDataset(:,2))) || ...
        any(cellfun(@ndims,modelDataset(:,2))-2) || ...
        any(cellfun('size',modelDataset(:,2),1)-1)
    errId = ['MAPS:',mfilename,':BadColumn2'];
    generate_MAPS_exception_and_add_as_cause(ValidationE,errId);
end
if ~isempty(ValidationE.cause)
    throw(ValidationE);
end

%% CHECK VALIDITY OF METADATA
% Use the MAPS unpack model and lookup model meatdata function to check
% that the metadata in the cell array model dataset is part of the model
% input (and is a unique set of metadata).
try
    unpack_model_metadata_and_lookup_index_numbers(...
        Model,modelMetadataStr,modelDataset(:,1));
catch ModelIndexE
    errId = ['MAPS:',mfilename,':BadMetadata'];
    BadMetadataE = generate_MAPS_exception_and_add_cause(...
        ModelIndexE,errId);
    ValidationE = addCause(ValidationE,BadMetadataE);
end

%% CHECK VALIDITY OF DATA
% Convert the numeric row vectors of data in the second column of the cell
% array to an equivalent matrix of data. Call a MAPS helper to validate
% that that dataset contains only non-inf, real numbers (it may contain
% NaN, used to signify periods where there is no data for that variable in
% MAPS cell array incomplete model datasets).
modelDatasetMat = convert_cell_array_of_vectors_to_matrix_equivalent(...
    modelDataset(:,2));
ValidationE = validate_incomplete_model_dataset(...
    modelDatasetMat,modelDataset(:,1),ValidationE);

%% THROW ANY EXCEPTIONS CAUGHT
% If the master exception contains any causes, throw it back out of this
% function.
if ~isempty(ValidationE.cause)
    throw(ValidationE);
end

end