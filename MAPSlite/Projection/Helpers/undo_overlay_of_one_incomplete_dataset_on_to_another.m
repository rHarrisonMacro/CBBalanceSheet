function [modelDataset1,modelMetadata1,modelDataset2,modelMetadata2] = ...
    undo_overlay_of_one_incomplete_dataset_on_to_another(...
    modelDataset,modelMetadata,data1Logicals,data2Logicals)
% This helper un-does an overlay of one incomplete dataset on to another.
% It separates a combined incomplete dataset into ist constituent parts.
%
% INPUTS:
%   -> modelDataset: combined matrix of data   
%   -> modelMetadata: cell string array of model metadata or index numbers 
%      for that data
%   -> data1Logicals: matrix of logicals describing which elements belong 
%      to the first dataset   
%   -> data2Logicals: matrix of logicals describing which elements belong 
%      to the second dataset 
%
% OUTPUTS:  
%   -> modelDataset1: first dataset separated out  
%   -> modelMetadata1: cell string array of model metadata or index numbers 
%      for that data
%   -> modelDataset2: second dataset separated out   
%   -> modelMetadata2: cell string array of model metadata or index numbers 
%      for that data
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%
% DETAILS:  
%   -> This forecast macro helper separates two incomplete datasets that
%      have been combined into one for a particular variable type in the 
%      model (eg model variables).
%   -> The dataset is assumed to be potentially incomplete, meaning that
%      missing values may be represented by NaNs. In addition, the dataset
%      need not include all the variables of a given type (just those 
%      described in the metadata passed in).
%   -> This operation reverses the operation to combine two incomplete
%      datasets into one. See the documentation on that function for more
%      details.
%
% NOTES:   
%   -> See <> for a description of MAPS forecast macro helpers.
%   -> See also <> which reverses the operation in this function to
%      separate out the data input.
%
% This version: 07/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that four inputs were passed in and then use a helper function to
% check that the incomplete model datasets and metadata passed in are
% valid.
if nargin < 4
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~is_non_inf_real_two_dimensional_numeric_matrix(modelDataset)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~is_column_cell_string_array(modelMetadata) && ...
        ~is_finite_real_numeric_column_vector(modelMetadata)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId); 
elseif ~is_two_dimensional_logical_matrix(data1Logicals)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
elseif ~is_two_dimensional_logical_matrix(data2Logicals)
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);
end

%% CHECK DIMENSIONS OF DATA FOR CONSISTENCY
% Compute the dimensions of the dataset input. Throw an exception if the
% metadata and logical matrices input are not of the same dimenions.
[nData,S] = size(modelDataset);
if size(modelMetadata,1) ~= nData
    errId = ['MAPS:',mfilename,':MetadataDimsInconsistentWithData'];
    generate_and_throw_MAPS_exception(errId);    
elseif size(data1Logicals,1)~=nData || size(data1Logicals,2)~=S
    errId = ['MAPS:',mfilename,':Logical1DimsInconsistentWithData'];
    generate_and_throw_MAPS_exception(errId);
elseif size(data2Logicals,1)~=nData || size(data2Logicals,2)~=S
    errId = ['MAPS:',mfilename,':Logical2DimsInconsistentWithData'];
    generate_and_throw_MAPS_exception(errId);    
end

%% CREATE A MATRIX OF LOGICALS
% Create a matrix of logicals for each of the individual datasets of the 
% same dimension as that dataset (which can be inferred by removing rows 
% that consistent entirely of zeros in the logicals input).
data1SepLogicals = data1Logicals;
data2SepLogicals = data2Logicals;
data1SepLogicals(~any(data1Logicals,2),:) = [];
data2SepLogicals(~any(data2Logicals,2),:) = [];

%% CREATE EMPTY MATRICES FOR THE SEPARATED DATASETS
% Create empty matrices (of NaNs) for the separated datasets by measuring 
% the dimensions of the logicals created above.
ndata1 = size(data1SepLogicals,1);
ndata2 = size(data2SepLogicals,1);
modelDataset1 = NaN*ones(ndata1,S);
modelDataset2 = NaN*ones(ndata2,S);

%% INSERT THE SEPARATED ELEMENTS INTO THE DATASETS
% Use the logicals computed above to insert the correct elements of the
% combined dataset into the separated out matrices. Compute the metadata
% for each dataset using the index numbers of the combined set with
% any rows of the input logicals that contained non-zeros. 
modelDataset1(data1SepLogicals) = modelDataset(data1Logicals);
modelDataset2(data2SepLogicals) = modelDataset(data2Logicals);
modelMetadata1 = modelMetadata(any(data1Logicals,2));
modelMetadata2 = modelMetadata(any(data2Logicals,2));

end