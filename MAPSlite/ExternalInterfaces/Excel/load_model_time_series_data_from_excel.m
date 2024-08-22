function varargout = load_model_time_series_data_from_excel(...
    Model,modelVarIds,excelFileName,sheetName,Options)
% This helper loads in model time series data from an Excel file.
% In particular, it loads in a complete set of time series for one or more
% variable types residing in a MAPS model structure.
%
% INPUTS:
%   -> Model: MAPS model object
%   -> modelVarIds: one or more identifiers to a MAPS variable type
%   -> excelFileName: name of the Excel file
%   -> sheetName (optional): name of the sheet in Excel file (can be '')
%   -> Options (optional): structure of options
%       - isComplete: true if data should not contain (NaN) holes
%       - isStraightEdge: true if all data should have same time dimension
%       - containsHeaders: true if headers are expected to be part of data
%       - isCompleteSetOfHeaders: true if all data should have headers
%       - containsDates: true if dates are expected to be part of dataset
%       - isRealDates: true if real dates are used (eg 31/12/2012) with the
%         alternative format being index numbers (i.e. 1 to T)
%
% OUTPUTS:
%   -> varagout: nVar*T matrix for each model identifier passed in
%
% DETAILS:
%   -> This function loads in time series data for one or more sets of 
%      variables in a MAPS model structure from an Excel spreadsheet.
%   -> The Excel spreadsheet can be formatted in a variety of different
%      ways as controlled by the setting of the options above.
%   -> The resulting dataset is formatted consistenly with time series data
%      throughout MAPS.
%
% NOTES:
%   -> This is part of a family of functions in MAPS which load in and save
%      out data.
%
% This version: 05/11/2012
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~ischar(modelVarIds) && ...
        ~is_column_cell_string_array(modelVarIds) && ...
        ~is_row_cell_string_array(modelVarIds)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);    
end

%% MANAGE OPTIONAL INPUTS
if nargin < 5
    Options = struct;
end
if nargin < 4
    sheetName = '';
end

%% DEAL WITH ALTERNATIVE SHAPE OF MODEL VAR IDS
if ischar(modelVarIds)
    modelVarIds = {modelVarIds};
elseif size(modelVarIds,2) == 1
    modelVarIds = modelVarIds';
end

%% CHECK NUMBER OF OUTPUT ARGUMENTS
nModelVarIds = size(modelVarIds,2);
if nargout > nModelVarIds+1
    errId = ['MAPS:',mfilename,':TooManyOutputArgs'];
    generate_and_throw_MAPS_exception(errId);     
end

%% GET MODEL VARAIBLE MENONICS
% This unpacks the mnemonics for each of the variable types requested in 
% modelVarIds input and then stacks them together into a single column cell
% string array containing all the mnemonics.
modelVarMnems = cell(1,nModelVarIds);
for iId = 1:nModelVarIds
    modelVarMnems{iId} = unpack_model(Model,modelVarIds(iId));
end
allModelVarMnems = vertcat(modelVarMnems{:});

%% LOAD THE TIME SERIES
% This helper function will load in all of the time series data and then
% extract the variables specified to return a single matrix of time series
% data ordered consistently with the variable mnemonics.
[dataMat,~,dataDateNums] = load_specified_time_series_data_from_excel(...
    allModelVarMnems,excelFileName,sheetName,Options);

%% EXTRACT THE SERIES IN INDIVIDUAL MATRICES
% Partition the data matrix returned above into separate matrices; one for
% each of the variable types. It combines information about the number of
% each variable with the knowledge that the matrix of time series data
% returned above will be ordered by variable and by each type.
varargout = cell(1,nModelVarIds+1);
nSeriesExtracted = 0;
for iId = 1:nModelVarIds
    niIdSeries = size(modelVarMnems{iId},1);
    indexStart = nSeriesExtracted+1;
    nSeriesExtracted = nSeriesExtracted+niIdSeries;
    varargout{iId} = dataMat(indexStart:nSeriesExtracted,:);
end
varargout{nModelVarIds+1} = dataDateNums;

end