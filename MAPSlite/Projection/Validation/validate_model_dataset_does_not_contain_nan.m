function validate_model_dataset_does_not_contain_nan(...
    modelDataset,modelMetadata)
% This helper validates that a model dataset contains only non-nan values.
% If it does find nan values it throws an exception detailing metadata 
% (mnemonics or names) for the series that contain nan (NaN) values.
%
% INPUTS:   
%   -> modelDataset: matrix of data for a particular model variable type   
%   -> modelMetadata: cell string array of model metadata
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause
%
% DETAILS:  
%   -> This forecast macro helper validates that a dataset for a particular 
%      variable type in the model (eg model variables) does not contain nan 
%      values. 
%   -> This ia a useful check to make if, for example, any of the model 
%      operations might include a '0/0' (eg if a raw observable data series
%      is 0 and a scalar time-varying trend included in the data
%      tranformation for that series is also 0).
%   -> If all the data passed in is non-nan (and the number / shape of the 
%      inputs is as expected), then this function returns no output.
%   -> If some of the data is nan, then this function constructs an
%      exception.
%   -> For each row of the dataset that contains a nan value, this helper 
%      adds an exception as cause with the metadata associated with that 
%      series in the message. 
%   -> The metadata is taken as the element corresponding to the row in the 
%      cell string array of metadata input.
%
% NOTES:   
%   -> See <> for a description of data errors exception handling in MAPS..
%
% This version: 22/02/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of the inputs is as expected. All inputs
% are compulsory.
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isnumeric(modelDataset) || ndims(modelDataset)~=2
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~iscellstr(modelMetadata)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
end

%% CHECK CONSISTENCY OF DATASET WITH METADATA
% Compute the number of series in the dataset (using the MAPS convention
% that model datasets are nSeries*T). Check that that is consistent with
% the metadata input. If not, throw an exception.
nSeries = size(modelDataset,1);
if size(modelMetadata,1) ~= nSeries
    errId = ['MAPS:',mfilename,':MetadataDimsInconsistentWithData'];
    generate_and_throw_MAPS_exception(errId);
end

%% SEARCH FOR NAN VALUES
% If any nan components are found (isnan is an element-by-element test), 
% create an exception to add individual causes to. Unpack the metadata
% from the model and check its dimension against the dataset. Throw an 
% exception if they're not compatible. Otherwise, add the metadata for each
% of the series found as causes to the exception.
if any(any(isnan(modelDataset)))
    masterErrId = ['MAPS:',mfilename,':BadDataset'];
    BadDatasetE = generate_MAPS_exception(masterErrId);
    errId = [masterErrId,':BadSeries'];
    indBadSeries = find(any(isnan(modelDataset),2));
    nBadSeries = size(indBadSeries,1);
    for iSeries = 1:nBadSeries
        BadDatasetE = generate_MAPS_exception_and_add_as_cause(...
            BadDatasetE,errId,modelMetadata(indBadSeries(iSeries)));
    end
    throw(BadDatasetE);
end

end