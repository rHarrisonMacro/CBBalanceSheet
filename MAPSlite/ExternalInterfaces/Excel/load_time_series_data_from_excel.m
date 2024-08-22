function [dataMat,dataHeaders,dateNums,dateStrs] = ...
    load_time_series_data_from_excel(excelFileName,sheetName,Options)
% This helper loads in time series data from a sheet in an Excel file.
%
% INPUTS:
%   -> excelFileName: name of the Excel file
%   -> sheetName (optional): name of the sheet in Excel file (can be '')
%   -> Options (optional): structure of options for the time series
%       - isComplete: true if data should not contain (NaN) holes
%       - isStraightEdge: true if all data should have same time dimension
%       - containsHeaders: true if headers are expected to be part of data
%       - isCompleteSetOfHeaders: true if all data should have headers
%       - isUniqueSetOfHeaders: true if all headers should be unique
%       - containsDates: true if dates are expected to be part of dataset
%       - isRealDates: true if real dates are used (eg 31/12/2012) with the
%         alternative format being index numbers (i.e. 1 to T)
%       - isStrDates: true if dates are in string format (eg 'dd/mmm/yy')
%         with tha laternative format being a numeric format
%       - strDateFormat: format for string dates (either in Excel file or
%         to be created if the dates in the file are numeric)
%
% OUTPUTS:
%   -> dataMat: nVars*T matrix of time series data
%   -> dataHeaders: nVars*1 column cell string array of headers
%   -> dateNums: 1*T row array of date numbers
%   -> dateStrs: 1*T row array of date strings (will be empty if
%      isRealDates is false)
%
% DETAILS:
%   -> This function loads in time series data residing in the specified 
%      worksheet of the specified file.
%   -> It does so in two steps using two helper functions: first, it loads
%      the raw data in; second, it cleans (and checks) it. See the content
%      of each function for details of each step.
%
% NOTES:
%   -> This is part of a family of functions in MAPS which load in and save
%      out data.
%
% This version: 01/11/2012
% Author(s): Matt Waldron

%% CHECK INPUTS
% Note that checking of the input types is left to the inner function
% called.
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
end

%% MANAGE OPTIONAL INPUTS
% This allows calls to the functions below to be independent of the
% optional set of input arguments passed into this function.
if nargin < 3
    Options = struct;
end
if nargin < 2
    sheetName = '';
end

%% LOAD RAW DATA
rawDataCell = load_raw_data_from_excel(excelFileName,sheetName);

%% CONVERT RAW LOADED DATA TO TIME SERIES
[dataMat,dataHeaders,dateNums,dateStrs] = ...
    extract_time_series_from_loaded_cell_array_of_data(...
    rawDataCell,Options);

end