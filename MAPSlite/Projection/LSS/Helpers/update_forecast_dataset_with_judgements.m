function fDataNew = update_forecast_dataset_with_judgements(...
    modelMetadata,fData,jCell)
% This helper updates a forecast dataset for judgements.
% It takes the values in the judgements input and overlays them on to the
% forecast dataset input.
%
% INPUTS:   
%   -> modelMetadata: cell string array of model metadata
%   -> fData: numeric forecast dataset
%   -> jCell: two-column cell array of judgement metadata and values
%
% OUTPUTS:  
%   -> fDataNew: numeric forecast dataset updated with judgements
%
% CALLS:    
%   -> generate_and_throw_MAPS_exception
%   -> lookup_model_index_numbers
%   -> convert_cell_array_of_vectors_to_matrix_equivalent
%
% DETAILS:  
%   -> This helper updates a forecast dataset for judgements.
%   -> It is useful in applications where judgements made can be applied 
%      directly to a particular forecast dataset (eg judgements made 
%      directly to shocks in linear state space models and directly to 
%      residuals in non-linear backward-looking models). 
%   -> The operation is completed in three steps. First, the judgements are
%      associated with variables in the dataset through an index lookup of
%      the judgement metadata in the model metadata. Second, the individual 
%      vectors of judgement data are combined into a single matrix (where 
%      missing periods are filled in with NaNs). Finally, the judgements
%      are overlaid on to the input dataset.
%   -> Other than input and input consistency errors, this function will
%      throw an exception if the horizon of the judgements exceeds that of
%      the dataset input.
%           
% NOTES:    
%   -> Please see <> for a description of MAPS module helpers.
%   -> See also <> for a description of imposing judgement in MAPS.
%           
% This version: 04/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of the inputs is as expected. All inputs
% are compulsory.
if nargin < 3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId);
elseif ~iscellstr(modelMetadata) || ndims(modelMetadata)~=2 || ...
        size(modelMetadata,2)~=1
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);  
elseif ~(isnumeric(fData)&&ndims(fData)==2)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif ~iscell(jCell) || ndims(jCell)~=2 || size(jCell,2)~=2
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
end

%% CHECK CONSISTENCY OF MEATDATA & DATASET
% Check that the dataset input has dimensions consistent with the model
% metadata input. If not, throw an exception.
nfData = size(fData,1);
if size(modelMetadata,1) ~= nfData
    errId = ['MAPS:',mfilename,':DatasetIncompatibleWithMetadata'];
    generate_and_throw_MAPS_exception(errId);
end

%% LOOKUP INDEX NUMBERS OF JUDGEMENTS IN MODEL
% Use a MAPS utility function to lookup the index numbers of the judgements
% in the model metadata passed in.
try
    jInds = lookup_model_index_numbers(modelMetadata,jCell(:,1));
catch MetadataLookupE
    errId = ['MAPS:',mfilename,':IndexLookupFailure'];
    generate_MAPS_exception_add_cause_and_throw(MetadataLookupE,errId);
end

%% CONVERT CELL ARRAY TO MATRIX EQUIVALENT
% Use another MAPS utility function to convert the cell array of numeric
% row vectors to a matrix equivalent.
try
    jMat = convert_cell_array_of_vectors_to_matrix_equivalent(jCell(:,2));
catch DataConversionE
    errId = ['MAPS:',mfilename,':CellDataConversionFailure'];
    generate_MAPS_exception_add_cause_and_throw(DataConversionE,errId);
end

%% CHECK JUDGEMENT TIME HORIZON IS COMPATIBLE WITH DATA
% Compute the horizon of the forecast dataset input and the judgement
% matrix computed above. If the judgement horizon is larger than the
% dataset horizon, throw an exception.
H = size(fData,2);
Hj = size(jMat,2);
if Hj > H
    errId = ['MAPS:',mfilename,':JudgementHorizonError'];
    generate_and_throw_MAPS_exception(errId);
end

%% SETUP OUTPUT
% Set the forecast dataset output equal to the dataset input.
fDataNew = fData;

%% UPDATE DATASET VALUES
% Determine the non-NaN components of the judgement matrix. Setup a logical
% matrix to be used to replace elements of the forecast dataset. Convert
% the elements in that matrix that correspond to non-NaNs to  
% Setup a logical matrix to be used below to replace values in fData.
% Replace those elements with the non-NaN values in the judgement matrix.
[sortedIndsForfData,sortedIndsForjMat]=sort(jInds);
jMat = jMat(sortedIndsForjMat,:);
nonNaNlogicals = ~isnan(jMat);
fDataToReplaceLogicals = false(nfData,H);
fDataToReplaceLogicals(sortedIndsForfData,1:Hj) = nonNaNlogicals;
fDataNew(fDataToReplaceLogicals) = jMat(nonNaNlogicals);

end