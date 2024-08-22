function [dataCell,colHeaders] = ...
    load_raw_data_from_excel_with_column_headers(excelFileName,sheetName)
% This helper loads in data from an Excel file with column headers.
%
% INPUTS:
%   -> excelFileName: name of the Excel file
%   -> sheetName (optional): name of the sheet in the Excel file
%
% OUTPUTS:
%   -> dataCell: cell array of loaded data (with headers separated out)
%   -> colHeaders: cell string array of column headers for that data
%
% DETAILS:
%   -> This function loads in all data residing in the specified worksheet
%      of the specified file.
%   -> It then extracts the first row of that data as column headers.
%   -> It allows for the possibility that some column headers might be
%      empty, but will error if non-empty components of the 1st row of the
%      spreadsheet do not consitute valid strings in MATLAB.
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

%% LOAD RAW DATA
if nargin > 1
    rawDataCell = load_raw_data_from_excel(excelFileName,sheetName);
else
    rawDataCell = load_raw_data_from_excel(excelFileName);
end

%% EXTRACT 1ST ROW FROM REST
rawColHeaders = rawDataCell(1,:);

%% CLEAN AND CHECK COLUMN HEADERS
colHeaders = clean_column_headers_loaded_from_excel(rawColHeaders);

%% EXTRACT REST OF THE DATA
% Extract all the data under the headers, checking that there is any!
if size(rawDataCell,1) < 2
    errId = ['MAPS:',mfilename,':NoDataUnderHeaders'];
    generate_and_throw_MAPS_exception(errId);
end
dataCell = rawDataCell(2:end,:);

end