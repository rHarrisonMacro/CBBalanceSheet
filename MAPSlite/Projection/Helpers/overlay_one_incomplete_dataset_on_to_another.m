function [modelDataset,modelMetadata,data1Logicals,data2Logicals] = ...
    overlay_one_incomplete_dataset_on_to_another(...
    modelDataset1,modelMetadata1,modelDataset2,modelMetadata2)
% This helper overlays one incomplete dataset on to another.
% It combines the two incomplete datasets (with the first taking
% precedence) using the metadata information passed in for each dataset.
%
% INPUTS:
%   -> modelDataset1: matrix of data to overlay   
%   -> modelMetadata1: cell string array of model metadata or index numbers 
%      for that data
%   -> modelDataset2: matrix of data to be overlaid on   
%   -> modelMetadata2: cell string array of model metadata or index numbers 
%      for that data
%
% OUTPUTS:  
%   -> modelDataset: combined matrix of data   
%   -> modelMetadata: cell string array of model metadata for that data
%   -> data1Logicals: matrix of logicals describing which elements belong 
%      to the first dataset   
%   -> data2Logicals: matrix of logicals describing which elements belong 
%      to the second dataset 
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%
% DETAILS:  
%   -> This forecast macro helper combines two incomplete datasets into one
%      a particular variable type in the model (eg model variables).
%   -> The dataset is assumed to be potentially incomplete, meaning that
%      missing values may be represented by NaNs. In addition, the dataset
%      need not include all the variables of a given type (just those 
%      described in the metadata passed in).
%   -> If any of the metadata overlaps in the two datasets, then this
%      function will combine that data into one (with the data in the first
%      dataset taking precedence if it is non-NaN).
%   -> The operation is executed on a period-by-period basis so that the
%      resulting combined matrix can contain rows (time series of data for 
%      a particular variable from the model) that are the combination of 
%      data from both datasets input.
%
% NOTES:   
%   -> See <> for a description of MAPS forecast macro helpers.
%   -> See also <> which reverses the operation in this function to
%      separate out the combined set of data that is output to this 
%      function into the original inputs.
%
% This version: 07/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that four inputs were passed in and then check that the incomplete 
% model datasets and metadata passed in are valid.
if nargin < 4
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~is_non_inf_real_two_dimensional_numeric_matrix(modelDataset1)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~is_column_cell_string_array(modelMetadata1) && ...
        ~is_finite_real_numeric_column_vector(modelMetadata1)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId); 
elseif ~is_non_inf_real_two_dimensional_numeric_matrix(modelDataset2)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);     
elseif ~is_column_cell_string_array(modelMetadata2) && ...
        ~is_finite_real_numeric_column_vector(modelMetadata2)
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);        
end
   
%% CHECK DIMENSIONS OF DATA FOR CONSISTENCY
% Compute the number of periods in the first incomplete model dataset
% input. Throw an exception if the seond dataset is not of the same
% dimension. Compute the number of variables in the first and second
% datasets. Throw an exception if the metadata for those datasets does not
% have consistent dimensions. 
[nData1,S] = size(modelDataset1);
nData2 = size(modelDataset2,1);
if size(modelDataset2,2) ~= S
    errId = ['MAPS:',mfilename,':Dataset1DimsInconsistentWithDataset2'];
    generate_and_throw_MAPS_exception(errId);
elseif size(modelMetadata1,1) ~= nData1
    errId = ['MAPS:',mfilename,':Metadata1DimsInconsistentWithData1'];
    generate_and_throw_MAPS_exception(errId);    
elseif size(modelMetadata2,1) ~= nData2
    errId = ['MAPS:',mfilename,':Metadata2DimsInconsistentWithData2'];
    generate_and_throw_MAPS_exception(errId);    
end

%% DETERMINE UNION OF THE METADATA
% Find the unique set of model meatdata passed in and the index numbers of
% that unique set.
[modelMetadata,uniqueMetadataInds] = unique(...
    [modelMetadata1;modelMetadata2]);

%% SETUP BASELINE MATRIX
% Using the unique indices produced above, setup a baseline matrix. If
% there were no overlap in the input vectors of model index numbers, then 
% this matrix is the same as the eventual output. If there was overlap,
% then the unique function assigns the indices to the second of the two
% repeat indices.
modelDataset12 = [modelDataset1;modelDataset2];
modelDataset = modelDataset12(uniqueMetadataInds,:);

%% SET THE LOGICAL OUTPUT VALUES
% Determine the dimensionality of all the inputs. Setup a matrix of
% logicals which show which elements in the combined matrix belong to each
% of the inputs. Set all of the indexes of the indices computed above that
% that fall below or equal to the number of rows in the first dataset and 
% all of the indexes above that to the second dataset.
nData = size(modelDataset,1);
data1Logicals = false(nData,S);
data2Logicals = false(nData,S);
data1Logicals(uniqueMetadataInds<=nData1,:) = true;
data2Logicals(uniqueMetadataInds-nData1>0,:) = true;

%% OVERWRITE THE INTERSECTION WITH NON-NAN VALUES FROM THE FIRST DATASET
% Because overlapping data was taken from the second dataset input, it is
% necessary to overwrite these numbers with values from the first matrix. 
% Find the intersecting indices and loop through them assigning each non-
% NaN number found in the first dataset as belonging to that dataset.
[intersectMetadata,indData1Reps,~] = intersect(...
    modelMetadata1,modelMetadata2);
[~,indDataReps,~] = intersect(modelMetadata,intersectMetadata);
nDataReps = size(indDataReps,1);
for iRep = 1:nDataReps
    idataRepToAdd = ~isnan(modelDataset1(indData1Reps(iRep),:));
    modelDataset(indDataReps(iRep),idataRepToAdd) = modelDataset1(...
        indData1Reps(iRep),idataRepToAdd);
    data1Logicals(indDataReps(iRep),idataRepToAdd) = true;
    data2Logicals(indDataReps(iRep),idataRepToAdd) = false;
end

end