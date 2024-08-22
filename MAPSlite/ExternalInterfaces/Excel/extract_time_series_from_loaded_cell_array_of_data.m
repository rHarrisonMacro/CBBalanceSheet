function [dataMat,dataHeaders,dateNums,dateStrs] = ...
    extract_time_series_from_loaded_cell_array_of_data(...
    loadedDataCell,Options)
% This helper extracts time series data from a cell array of loaded data.
%
% INPUTS:
%   -> loadedDataCell: cell array of (scalar) data
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
%   -> This function extracts MAPS-compatible time series data from a cell
%      array of loaded data.
%   -> The data could have been loaded from anywhere, but is typically
%      likely to come from Excel or CSV files.
%   -> Time series in MAPS are of the dimension nVars*nPeriods and the
%      output of this function is consistent with that. By contrast, the
%      input data cell is assumed to be consistent with an nPeriods*nVars
%      format which is typical of how a spreadsheet would be setup.
%   -> The optional options input allows control over the expected format
%      of the data and associated error handling. Defaults for each option
%      are stored in this function - see below.
%
% NOTES:
%   -> This is part of a family of functions in MAPS which load in and save
%      out data.
%
% This version: 05/11/2012
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~is_two_dimensional_cell_array(loadedDataCell) || ...
        size(loadedDataCell,1)<2 || size(loadedDataCell,2)<2
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif nargin>1 && ~isstruct(Options)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
end

%% MANAGE OPTIONAL INPUT
% This allows the code below to flow through without knowledge of the
% number of input arguments passed in.
if nargin < 2
    Options = struct;
end

%% EXTRACT OPTIONS & CHECK THEY ARE VALID
% Note that it is not possible to specify that the dataset is incomplete
% and that it has a straight edge. A dataset with a ragged edge is a
% special form of incompleteness, so a straight edge implies completeness.
[isComplete,isStraightEdge,containsHeaders,isCompleteSetOfHeaders,...
    isUniqueSetOfHeaders,containsDates,isRealDates,isStrDates,...
    strDateFormat] = get_time_series_data_options(Options);

%% REMOVE TRAILING NANS FROM DATA CELL
% This is useful because in some instances even after data have been 
% seemingly deleted from a cell in Excel, a trace of it remains, which 
% MATLAB reads in as NaNs.
cleanedLoadedDataCell = remove_trailing_NaNs(...
    loadedDataCell,containsHeaders);

%% EXTRACT SEPARATE PARTS OF THE DATA CELL
% This splits the data cell into dates (if appropriate for the option), 
% time series data, and headers.
[dataCell,colHeaders,dates] = ...
    extract_components_of_loaded_data_cell(...
    cleanedLoadedDataCell,containsDates,containsHeaders);

%% CONVERT CELL ARRAY OF NUMERIC DATA TO MATRIX EQUIVALENT
dataMat = convert_cell_array_of_numerics_to_matrix(...
    dataCell,isComplete,isStraightEdge);

%% CREATE DATA HEADERS FROM COLUMN HEADERS
dataHeaders = clean_column_headers(colHeaders,containsHeaders,...
    isCompleteSetOfHeaders,isUniqueSetOfHeaders);

%% CREATE DATES OUTPUT
% If the dates are supposed to be "real" (eg 31/12/2012) then check that
% the MATLAB datenum function recognises them. If not convert the dates
% cell to a matrix.
[dateNums,dateStrs] = create_dates_output(...
    dates,containsDates,isRealDates,isStrDates,strDateFormat,dataMat);

end

%% FUNCTION TO CREATE DATES OUTPUT FROM 1ST COLUMN OF DATA CELL ARRAY
function [dateNums,dateStrs] = create_dates_output(...
    dates,containsDates,isRealDates,isStrDates,strDateFormat,dataMat)
% This helper converts the dates input into string & numeric output format.
%
% INPUTS:
%   -> dates: row cell array of dates loaded from 1st column of cell array
%   -> containsDates: true/false
%   -> isRealDates: true/false
%   -> isStrDates: true/false
%   -> strDateFormat: eg 'dd/mm/yyyy'
%   -> nVars*T matrix of time series data
%
% OUTPUTS:
%   -> dateNums: 1*T vector of numeric dates representation
%   -> dateStrs: 1*T cell array of string dates representation 

%% COMPUTE NUMBER OF TIME PERIODS IN TIME SERIES MATRIX
nTimePeriods = size(dataMat,2);

%% SET OUTPUT FOR THE CASE IN WHICH NO DATES WERE IN DATASET
% If the dataset does not contain dates setup an index as output for the
% numeric dates and set the string dates to empty.
if ~containsDates
    dateStrs = [];    
    dateNums = (1:nTimePeriods);
    return
end

%% SET OUTPUT ACCORDING TO WHETHER DATES WERE REAL OR OTHERWISE
% If the dataset contains real dates set the string dates output to the
% dates input and compute the date numerics using the MATLAB datenum 
% function. Otherwise, check that the date index in the datset is valid and
% set the date strings to empty. Note that in both cases the transpose of
% the dates loaded in is taken so that the dates match the MAPS convention
% of time periods being read across columns.
if isRealDates
    [dateNums,dateStrs] = convert_real_dates_format(...
        dates',isStrDates,strDateFormat);    
else
    [dateNums,dateStrs] = convert_index_dates_format(dates',nTimePeriods);
end

end

%% FUNCTION TO DEAL WITH REAL DATES FORMAT
function [dateNums,dateStrs] = convert_real_dates_format(...
    dates,isStrDates,strDateFormat)
% This helper converts real dates into the format required for output.
%
% INPUTS:
%   -> dates: row cell array of dates loaded from 1st column of cell array
%   -> isStrDates: true/false
%   -> strDateFormat: eg 'dd/mm/yyyy'
%
% OUTPUTS:
%   -> dateNums: 1*T vector of numeric dates representation
%   -> dateStrs: 1*T cell array of string dates representation 

%% CHECK FORMAT OF DATES
% Dates must either be numeric or string - they cannot be a mixture.
isAllNumeric = all(cellfun(@(x) isnumeric(x),dates));
isAllString = iscellstr(dates);
if ~isAllNumeric && ~isAllString
    errId = ['MAPS:',mfilename,':BadRealDateFormat'];
    generate_and_throw_MAPS_exception(errId);  
elseif isAllNumeric && isStrDates
    errId = ['MAPS:',mfilename,':BadRealStringDatesFormat'];
    generate_and_throw_MAPS_exception(errId); 
elseif isAllString && ~isStrDates
    errId = ['MAPS:',mfilename,':BadRealNumericDatesFormat'];
    generate_and_throw_MAPS_exception(errId); 
end

%% IF DATES ARE NUMERIC
% Convert the cell array of dates to a numeric vector and then create a
% date string allowing for the fact that the MATLAB base year is 0000,
% whereas Excel uses 1900 (with the -1 accounting for some sort of 
% difference in the count of leap years).
if isAllNumeric    
    dateNums = cell2mat(dates);
    MATLABbaseDateNum = datenum('01/01/1900','dd/mm/yyyy')-...
            datenum('01/01/0000','dd/mm/yyyy')-1;
    dateNums = dateNums+MATLABbaseDateNum;
    try      
        dateStrs = cellstr(datestr(dateNums,strDateFormat))';
    catch DatstrE
        errId = ['MAPS:',mfilename,':BadRealDateNumFormat'];
        generate_MAPS_exception_add_cause_and_throw(...
            DatstrE,errId,{strDateFormat});
    end
end

%% IF DATES ARE STRINGS
% Set the date strings equal to dates and then use the datenum function to
% convert to the numeric format.
if isAllString
    dateStrs = dates;
    try
        dateNums = datenum(dateStrs,strDateFormat)';
    catch DatenumE
        errId = ['MAPS:',mfilename,':BadRealDateStrFormat'];
        generate_MAPS_exception_add_cause_and_throw(...
            DatenumE,errId,{strDateFormat});
    end
end

end

%% FUNCTION TO DEAL WITH INDEX DATES FORMAT
function [dateNums,dateStrs] = convert_index_dates_format(...
    dates,nTimePeriods)
% This helper converts index dates into the format required for output.
%
% INPUTS:
%   -> dates: row cell array of dates loaded from 1st column of cell array
%   -> nTimePeriods: number of time periods
%
% OUTPUTS:
%   -> dateNums: 1*T vector of numeric dates representation
%   -> dateStrs: 1*T cell array of string dates representation 

%% SET DATE STRINGS TO EMPTY
dateStrs = [];

%% CHECK FORMAT OF DATES
% If index dates have been used then the dates cell array must be numeric.
if ~all(cellfun(@(x) isnumeric(x),dates))
    errId = ['MAPS:',mfilename,':BadIndexDateFormat'];
    generate_and_throw_MAPS_exception(errId);    
end  

%% CONVERT FORMAT OF DATES
% Convert the cell array to a matrix and then check that the indices run
% from 1 to nTimePeriods.
dateNums = cell2mat(dates);
if ~isequal(dateNums,(1:nTimePeriods))
    errId = ['MAPS:',mfilename,':BadIndexDates'];
    generate_and_throw_MAPS_exception(errId);
end

end

%% FUNCTION TO CREATE TIME SERIES DATA MATRIX
function dataMat = convert_cell_array_of_numerics_to_matrix(...
    dataCell,isComplete,isStraightEdge)
% This helper converts the cell array of times series data to a matrix.
%
% INPUTS:
%   -> dataCell: cell array of numeric scalars
%   -> isComplete: true/false
%   -> isStraightEdge: true/false
%
% OUTPUTS:
%   -> dataMat: nVars*T matrix of time series data

%% CONVERT CELL ARRAY OF NUMERIC DATA TO MATRIX EQUIVALENT
% Catch any exceptions and throw with an appropriate error message. Note
% that the data cell is transposed before being converted to a matrix. This
% ensures that it conforms with the MAPS convention that time periods are
% measured across columns.
try
    dataMat = cell2mat(dataCell');
catch DataConversionE
    errId = ['MAPS:',mfilename,':DataConversionFailure'];
    generate_MAPS_exception_add_cause_and_throw(...
        DataConversionE,errId); 
end

%% CHECK DATASET
% Throw an error if: data contains non-real or inf values; data contains
% NaNs (holes) part way through it (eg 1;NaN;2) & the dataset should be 
% complete; one or more of the data series contains trailing NaNs and the
% data should have a straight (non-ragged) edge.
[containsNonRealOrInfVals,containsHoles,containsNaNs] = ...
    check_time_series_dataset(dataMat);
if containsNonRealOrInfVals
    errId = ['MAPS:',mfilename,':BadData'];
    generate_and_throw_MAPS_exception(errId); 
elseif isComplete && containsHoles
    errId = ['MAPS:',mfilename,':IncompleteData'];
    generate_and_throw_MAPS_exception(errId);     
elseif isStraightEdge && containsNaNs
    errId = ['MAPS:',mfilename,':RaggedEdgeOnData'];
    generate_and_throw_MAPS_exception(errId);
end

end

%% FUNCTION TO CHECK PROPERTIES OF DATASET
function [containsNonRealOrInfVals,containsHoles,containsNaNs] = ...
    check_time_series_dataset(dataMat)
% This helper examines some of the properties of the time series data.
%
% INPUTS:
%   -> dataMat: nVars*T matrix of time series data
%
% OUTPUTS:
%   -> containsNonRealOrInfVals: true/false
%   -> containsHoles: true/false
%   -> containsNaNs: true/false

%% CHECK FOR NON-REAL OR INF VALUES
% Uses a MAPS helper function to check if any of the elements of the matrix
% are non-real or infinite.
containsNonRealOrInfVals = ~is_non_inf_real_two_dimensional_numeric_matrix(...
    dataMat);

%% CHECK FOR NANS
% Check to see if there are any missing values. These would not be allowed
% if the dataset should be complete and should have a straight edge.
dataIsNaN = isnan(dataMat);
containsNaNs = any(any(dataIsNaN));

%% CHECK FOR HOLES IN DATASET
% By assertaining which elements of the matrix are NaNs and then prefixing
% that with a col of "false" logicals, we can see if there are any holes in
% the dataset (making an allowance for the possibility that no data may 
% exist in a particular column). For example, consider cases for a single
% time series of data with 5 observations. The following are incomplete: 
% [1;NaN;3;4;5] & [NaN;NaN;3;4;5]. By contrast, the following are complete:
% [1;2;3;4;5], [1;2;3;NaN;NaN] & [NaN;NaN;NaN;NaN;NaN].
nDataSeries = size(dataMat,1);
dataIsNaNaugmented = [false(nDataSeries,1) dataIsNaN];
containsHoles = any(...
    any(diff(dataIsNaNaugmented,1,2)<0,2)&~all(dataIsNaN,2));

end

%% FUNCTION TO CLEAN COLUMN HEADERS
function dataHeaders = clean_column_headers(...
    colHeaders,containsHeaders,isCompleteSetOfHeaders,isUniqueSetOfHeaders)
% This helper cleans loaded column headers for output.
% It replaces any NaN column headers with empty strings and then checks 
% that the data headers constitute a row cell string array. 
%
% INPUTS:
%   -> colHeaders: row cell array of column headers
%   -> containsHeaders: true/false
%   -> isCompleteSetOfHeaders: true/false
%   -> isUniqueSetOfHeaders: true/false
%
% OUTPUTS:
%   -> dataHeaders: cell string array of data headers

%% SET DATA HEADERS TO EMPTY IF NON EXIST
if ~containsHeaders
    dataHeaders = [];
    return
end

%% REPLACE NAN COLUMN HEADERS WITH EMPTY STRING
% Clean the column headers and transpose them so that they match the MAPS
% convention of variables measured down rows and time periods measured
% across columns.
dataHeaders = clean_column_headers_loaded_from_excel(colHeaders)';

%% CHECK HEADERS
% Throw an error if some of them are empty when they should be complete or
% if some of them are repeated when they should be unique (allowing for the
% possibility that there may be multiple missing data headers of that 
% option is set to false).
containsEmptyStrs = any(strcmp('',dataHeaders));
if isCompleteSetOfHeaders && containsEmptyStrs
    errId = ['MAPS:',mfilename,':MissingDataHeaders'];
    generate_and_throw_MAPS_exception(errId);
end
nonUniqueHeaders = find_repeated_strings(dataHeaders(~containsEmptyStrs));
if isUniqueSetOfHeaders && ~isempty(nonUniqueHeaders)
    errId = ['MAPS:',mfilename,':NonUniqueDataHeaders'];
    generate_and_throw_MAPS_exception(errId);
end

end

%% FUNCTION TO SPLIT LOADED DATA CELL INTO CONSTITUENT PARTS
function [dataCell,colHeaders,dates] = ...
    extract_components_of_loaded_data_cell(...
    loadedDataCell,containsDates,containsHeaders)
% This helper splits a data cell into its constituent parts.
% It uses the expected MAPS format of externally loaded time series to
% split a loaded time series data cell into data, dates and data headers 
% (as apprporiate given the options invoked). 
%
% INPUTS:
%   -> loadedDataCell: loaded cell array of scalar data
%   -> containsDates: true/false
%   -> containsHeaders: true/false
%
% OUTPUTS:
%   -> dataCell: loaded cell array of scalars
%   -> colHeaders: cell array of column headers
%   -> dates: cell array of dates

%% DEFINE DATA ROW AND COLUMN STARTING INDEX
if containsDates
    dataStartCol = 2;
else
    dataStartCol = 1;
end
if containsHeaders
    dataStartRow = 2;
else
    dataStartRow = 1;
end
    
%% EXTRACT DATA AND DATES / HEADERS AS APPROPRIATE
dataCell = loadedDataCell(dataStartRow:end,dataStartCol:end);
if containsDates
    dates = loadedDataCell(dataStartRow:end,1);
else
    dates = [];
end
if containsHeaders
    colHeaders = loadedDataCell(1,dataStartCol:end);
else
    colHeaders = [];
end

end

%% FUNCTION TO REMOVE TRAILING NANS
function cleanedLoadedDataCell = remove_trailing_NaNs(...
    loadedDataCell,containsHeaders)
% This helper removes trailing NaNs from a loaded cell array of data.
% This is useful because sometimes the "memory" of data is left in Excel
% which MATLAB reads in as NaNs.
%
% INPUTS:
%   -> loadedDataCell: loaded cell array of scalars
%   -> containsHeaders: true/false
%
% OUTPUTS:
%   -> cleanedLoadedDataCell: cell array of scalars with trailing NaNs
%      removed

%% FIND ROWS COMPRISED ENTIRELY OF NANS
isNumeric = cellfun(@(x) isnumeric(x),loadedDataCell);
entireRowIsNumeric = all(isNumeric,2);
entireRowIsNaN = entireRowIsNumeric;
entireRowIsNaN(entireRowIsNumeric) = ...
    all(cellfun(@(x) isnan(x),loadedDataCell(entireRowIsNumeric,:)),2);

%% COUNT NUMBER OF ROWS ONE END OF DATASET WITH ALL NANS
nDataCellRows = size(loadedDataCell,1);
nTrailingNaNs = 0;
for iRow = nDataCellRows:-1:1
    if entireRowIsNaN(iRow)
        nTrailingNaNs = nTrailingNaNs+1;
    else
        break
    end
end

%% THROW ERROR IF IS NO DATA LEFT
if (containsHeaders&&(nTrailingNaNs>=nDataCellRows-1)) || ...
        (~containsHeaders&&(nTrailingNaNs==nDataCellRows))
    errId = ['MAPS:',mfilename,':EmptyTimeSeriesData'];
    generate_and_throw_MAPS_exception(errId);
end

%% REMOVE TRAILING NANS
cleanedLoadedDataCell = loadedDataCell(1:nDataCellRows-nTrailingNaNs,:);

end

%% FUNCTION TO DEFINE DEFAULT OPTIONS FOR LOADING OF TIME SERIES DATA
function [isComplete,isStraightEdge,containsHeaders,...
    isCompleteSetOfHeaders,isUniqueSetOfHeaders,containsDates,...
    isRealDates,isStrDates,strDateFormat] = ...
    get_time_series_data_options(UserOptions)
% This helper extracts time series data from a cell array of loaded data.
%
% INPUTS:
%   -> UserOptions: structure of user options for the time series (can be an
%      empty structure)
%
% OUTPUTS:
%   -> isComplete: true/false
%   -> isStraightEdge: true/false
%   -> containsHeaders: true/false
%   -> isCompleteSetOfHeaders: true/false
%   -> isUniqueSetOfHeaders: true/false
%   -> containsDates: true/false
%   -> isRealDates: true/false
%   -> isStrDates: true/false
%   -> strDateFormat: eg 'dd/mm/yyyy'

%% DEFINE DEFAULT OPTIONS
DefaultOptions.isComplete = true;
DefaultOptions.isStraightEdge = true;
DefaultOptions.containsHeaders = true;
DefaultOptions.isCompleteSetOfHeaders = true;
DefaultOptions.isUniqueSetOfHeaders = true;
DefaultOptions.containsDates = true;
DefaultOptions.isRealDates = true;
DefaultOptions.isStrDates = true;
DefaultOptions.strDateFormat = 'dd/mm/yyyy';

%% OVERLAY INPUT OPTIONS ON TO DEFAULTS
Options = overlay_default_structure(UserOptions,DefaultOptions);

%% UNPACK OPTIONS
isComplete = Options.isComplete;
isStraightEdge = Options.isStraightEdge;
containsHeaders = Options.containsHeaders;
isCompleteSetOfHeaders = Options.isCompleteSetOfHeaders;
isUniqueSetOfHeaders = Options.isUniqueSetOfHeaders;
containsDates = Options.containsDates;
isRealDates = Options.isRealDates;
isStrDates = Options.isStrDates;
strDateFormat = Options.strDateFormat;

%% CHECK OPTIONS DATA TYPE
if ~is_logical_scalar(isComplete) || ...
        ~is_logical_scalar(isStraightEdge) || ...
        ~is_logical_scalar(containsHeaders) || ...
        ~is_logical_scalar(isCompleteSetOfHeaders) || ...
        ~is_logical_scalar(containsDates) || ...
        ~is_logical_scalar(isRealDates) || ...
        ~is_logical_scalar(isUniqueSetOfHeaders) || ...
        ~is_logical_scalar(isStrDates) || ...
        ~ischar(strDateFormat)
    errId = ['MAPS:',mfilename,':BadOption'];
    generate_and_throw_MAPS_exception(errId);
end

%% CHECK COMBINATIONS OF OPTIONS
if isStraightEdge && ~isComplete
    errId = ['MAPS:',mfilename,':BadOptionCombination'];
    generate_and_throw_MAPS_exception(errId);
end   

end