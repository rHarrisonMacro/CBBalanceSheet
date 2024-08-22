function [dataMat,dataHeaders,dateNums,dateStrs] = ...
    load_specified_time_series_data_from_excel(...
    dataHeaders,excelFileName,sheetName,Options)
% This helper loads in specified time series data from an Excel file.
% It allows a subset of time series data from a sheet in an Excel file to
% be loaded in. It also re-orders that data to be consistent with the
% ordering of the variables input.
%
% INPUTS:
%   -> dataHeaders: string or cell string vector of identifiers
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
%   -> dataMat: nVars*T matrix of time series data
%   -> dataHeaders: nVars*1 column cell string array of headers
%   -> dateNums: 1*T row array of date numbers
%   -> dateStrs: 1*T row array of date strings (will be empty if
%      isRealDates is false)
%
% DETAILS:
%   -> This function loads in specified time series data residing in the 
%      specified worksheet of the specified file.
%   -> It loads in the time series data from the spreadsheet and then users
%      a helper function to pick out the particular series required in the
%      order desired.
%   -> The data output is consistent with the conventions used in MAPS for
%      time series data. It is an nVars*T matrix. This function will also
%      optionally output the headers/identifiers/labels that go with that
%      data (as requested on input) as an nVars*1 cell string array, and a
%      1*T vector of date numbers and date strings as applicable. See the
%      content of the functions called below for more details.
%
% NOTES:
%   -> This is part of a family of functions in MAPS which load in and save
%      out data.
%
% This version: 05/11/2012
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~ischar(dataHeaders) && ...
        ~is_column_cell_string_array(dataHeaders) && ...
        ~is_row_cell_string_array(dataHeaders)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
end

%% MANAGE OPTIONAL INPUTS
if nargin < 4
    Options = struct;
end
if nargin < 3
    sheetName = '';
end

%% LOAD DATA FROM EXCEL
[allDataMat,allDataHeaders,dateNums,dateStrs] = ...
    load_time_series_data_from_excel(excelFileName,sheetName,Options);

%% CHECK VALIDITY OF DATA HEADERS
% Throw an error if the data headers output is empty (which would be the
% case if the spreadsheet did not contain headers and 'containsHeaders' was
% set to false in the options), or if the headers are non-unique (because
% it might not be possible to identify the variable of interest).
if isempty(allDataHeaders)
    errId = ['MAPS:',mfilename,':MissingHeaders'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isempty(find_repeated_strings(allDataHeaders))
    errId = ['MAPS:',mfilename,':NonUniqueHeaders'];
    generate_and_throw_MAPS_exception(errId);
end

%% EXTRACT SPECIFIED TIME SERIES
dataMat = extract_specified_time_series_from_dataset(...
    allDataMat,allDataHeaders,dataHeaders);

end