function rawDataCell = load_raw_data_from_excel(excelFileName,sheetName)
% This helper loads in all data from a worksheet in an Excel file.
%
% INPUTS:
%   -> excelFileName: name of the Excel file
%   -> sheetName (optional): name of the sheet in the Excel file
%
% OUTPUTS:
%   -> rawDataCell: cell array of loaded data
%
% DETAILS:
%   -> This function loads in all data residing in the specified worksheet
%      of the specified file.
%   -> The sheet name input is optional. If not provided, the load function
%      will attempt to load data from the 1st sheet in the file.
%
% NOTES:
%   -> This is part of a family of functions in MAPS which load in and save
%      out data.
%   -> If the excel file is not .xlsx or if the name of that file passed in
%      as input does not include a '.xlsx' extension than this function
%      will throw an error.
%
% This version: 01/11/2012
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~ischar(excelFileName)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif nargin>1 && ~ischar(sheetName)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
end

%% MANAGE OPTIONAL INPUT
% If sheet name was not provided, then set it to an empty string which
% allows the rest of the code below to flow through as if a sheet name had
% been input.
if nargin < 2
    sheetName = '';
end

%% CHECK FILE EXISTS ON SEACRH PATH WITH XLSX EXTENSION
try
    check_file_exists(excelFileName,'xlsx');
catch FileExistenceE
    errId = ['MAPS:',mfilename,':InvalidExcelDataFile'];
    generate_MAPS_exception_add_cause_and_throw(...
        FileExistenceE,errId,{excelFileName});
end

%% READ DATA
try
    [~,~,rawDataCell] = xlsread(excelFileName,sheetName);
catch ReadE
    errId = ['MAPS:',mfilename,':UnableToReadData'];
    generate_MAPS_exception_add_cause_and_throw(...
        ReadE,errId,{sheetName excelFileName});
end

%% THROW ERROR IF NOTHING WAS LOADED
% If the sheet was empty, then MATLAB reads in a scalar with a NaN value.
if isscalar(rawDataCell) && isnan(rawDataCell)
    errId = ['MAPS:',mfilename,':SheetIsEmpty'];
    generate_MAPS_exception_add_cause_and_throw(...
        errId,{excelFileName sheetName});
end

end