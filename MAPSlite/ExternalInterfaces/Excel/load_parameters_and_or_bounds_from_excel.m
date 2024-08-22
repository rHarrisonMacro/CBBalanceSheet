function LoadedParameters = load_parameters_and_or_bounds_from_excel(...
    excelFileName,sheetName)
% This helper loads a set of parameter mnemonic-value pairs from Excel.
% Alternatively, or in addition, it will also load in 
%
% INPUTS:
%   -> excelFileName: name of the Excel file
%   -> sheetName (optional): sheet name in the Excel file
%
% OUTPUTS:
%   -> LoadedParameters: structure of loaded parameter info:
%       - mnemonics: column cell string array of parameter mnemonics
%       - values: vector of parameter values
%       - lowerBounds: vector of lower bounds
%       - upperBounds: vector of upper bounds
%
% DETAILS:
%   -> This helper loads in mnemonic-value pairs and/or bounds from Excel.
%   -> It is designed to load in scalar parameters (rather than time series
%      - there are other functions for that).
%   -> It will load in all the parameters in the specified file and sheet,
%      as well as any values and/or bounds that are in that sheet.
%   -> It is useful as part of MAPS estimation functionality to ge
%      initialisations for optimisation.
%
% NOTES:
%   -> This function is part of a family of functions that load data from
%      Excel.
%
% This version: 23/11/2012
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
end

%% LOAD DATA FROM EXCEL
% The precise action depends on whether a sheet name was passed in.
if nargin < 2
    [dataCell,colHeaders] = ...
        load_raw_data_from_excel_with_column_headers(excelFileName);
else
    [dataCell,colHeaders] = ...
        load_raw_data_from_excel_with_column_headers(...
        excelFileName,sheetName);
end

%% CHECK DIMENSIONS OF DATA
% The sheet should contain two to four columns, one for mnemonics and one
% to three of parameter values and or bounds.
if size(dataCell,2)>4 || size(dataCell,2)<2
    errId = ['MAPS:',mfilename,':BadParamExcelSheetDims'];
    generate_and_throw_MAPS_exception(errId);
end

%% GET PERMITTED HEADERS CONFIGURATION
PermittedHeaders = get_permitted_parameter_sheet_headers();

%% CHECK & UNPACK COMPULSORY MNEM0NICS COLUMN
% The mnemonics must appear in the first column.
permittedMnemonicHeaders = PermittedHeaders.mnemonics;
if ~any(strcmpi(colHeaders{1},permittedMnemonicHeaders))
    errId = ['MAPS:',mfilename,':BadMnemonicsColHeader'];
    generate_and_throw_MAPS_exception(errId,{permittedMnemonicHeaders});
elseif ~is_column_cell_string_array(dataCell(:,1))
    errId = ['MAPS:',mfilename,':BadMnemsCol'];
    generate_and_throw_MAPS_exception(errId);
elseif size(unique(dataCell(:,1)),1)~=size(dataCell(:,1),1)
    errId = ['MAPS:',mfilename,':RepeatedMnems'];
    generate_and_throw_MAPS_exception(errId);
end
LoadedParameters.mnemonics = dataCell(:,1);

%% CHECK & UNPACK PARAMETER VECTOR COLUMNS
% For each of the vectors in the cell array loaded above, check whether
% they are headed values, lower bounds or upper bounds. If yes, unpack them
% for output. Handle possible errors, including that the parameter vectors
% are not numerics, the header is not recognised or that there is more than
% one header of the same type.
permittedParamVecTypes = {'values','lowerBounds','upperBounds'};
HeaderAlreadyFound = setup_header_existence_structure(...
    permittedParamVecTypes);
nPermittedTypes = size(permittedParamVecTypes,2);
nParamVecs = size(dataCell,2)-1;
for iParamVec = 1:nParamVecs
    iParamVecHeaderFound = false;
    iParamVecHeader = colHeaders{iParamVec+1};
    iParamVecCell = dataCell(:,iParamVec+1);
    for iType = 1:nPermittedTypes
        if ~iParamVecHeaderFound
            iParamVecTypeToSearchFor = permittedParamVecTypes{iType};
            [LoadedParameters,HeaderAlreadyFound,iParamVecHeaderFound] =...
                search_for_and_extract_parameter_vector(...
                LoadedParameters,HeaderAlreadyFound,PermittedHeaders,...
                iParamVecHeader,iParamVecCell,iParamVecTypeToSearchFor);
        end
    end
    if ~iParamVecHeaderFound
        errId = ['MAPS:',mfilename,':UnrecognisedHeader'];
        generate_and_throw_MAPS_exception(errId,{iParamVecHeader});
    end
end

end

%% FUNCTION TO SEARCH FOR AND UNPACK VALUES AND BOUNDS
function [LoadedParameters,HeaderAlreadyFound,paramVecTypeFound] = ...
    search_for_and_extract_parameter_vector(...
    LoadedParameters,HeaderAlreadyFound,PermittedHeaders,...
    colHeader,dataCol,paramVecTypeToSearchFor)
% This function searches for and extracts either param values or bounds.
% It checks that the specific type of vector - either parameter values,
% lower bounds or upper bounds - is contained in the input data column and
% then extracts it in numeric vector format as output to this function.
%
% INPUTS:
%   -> LoadedParameters: structure of outputs to this function
%       - mnemonics: column cell string array of parameter mnemonics
%       - values: vector of parameter values
%       - lowerBounds: vector of lower bounds
%       - upperBounds: vector of upper bounds
%   -> HeaderAlreadyFound: structure of logicals
%       - values: true/false
%       - lowerBounds: true/false
%       - upperBounds: true/false
%   -> PermittedHeaders: structure of permitted headers for:
%       - mnemonics: cell string array of permitted headers
%       - values: cell string array of permitted headers
%       - lowerBounds: cell string array of permitted headers
%       - upperBounds: cell string array of permitted headers
%   -> colHeader: string header above the column of data being considered
%   -> dataCol: column of data loaded from Excel
%   -> paramVecTypeToSearchFor: type of parameter to look for (either
%      values, lower bounds or upper bounds)
%
% OUTPUTS:
%   -> LoadedParameters: updated structure of outputs to this function
%       - mnemonics: column cell string array of parameter mnemonics
%       - values: vector of parameter values
%       - lowerBounds: vector of lower bounds
%       - upperBounds: vector of upper bounds
%   -> HeaderAlreadyFound: updated structure of logicals
%       - values: true/false
%       - lowerBounds: true/false
%       - upperBounds: true/false
%   -> paramVecTypeFound: true/false

%% IF PARAMETER TYPE COULD NOT BE FOUND, EXIT THIS FUNCTION
if ~any(strcmpi(colHeader,PermittedHeaders.(paramVecTypeToSearchFor)))
    paramVecTypeFound = false;
    return
end

%% IF PARAMETER TYPE HAS ALREADY BEEN FOUND, THROW AN ERROR
if HeaderAlreadyFound.(paramVecTypeToSearchFor)
    errId = ['MAPS:',mfilename,':RepeatedHeader'];
    generate_and_throw_MAPS_exception(errId,{colHeader});
end

%% SET LOGICALS FOR PARAMETER VECTOR TYPE
paramVecTypeFound = true;
HeaderAlreadyFound.(paramVecTypeToSearchFor) = true;

%% ATTEMPT TO CONVERT CELL ARRAY TO NUMERIC VECTOR
try
    LoadedParameters.(paramVecTypeToSearchFor) = cell2mat(dataCol);
catch Cell2MatE
    errId = ['MAPS:',mfilename,':BadValsCol'];
    generate_MAPS_exception_add_cause_and_throw(...
        Cell2MatE,errId,{colHeader});
end

end

%% FUNCTION TO GET PERMITTED HEADERS FOR SHEET COLUMNS
function PermittedHeaders = get_permitted_parameter_sheet_headers()
% This configuration defines valid column headers for a param excel sheet.
% Note that the check is not case sensitive, so all capitalised versions of
% the permitted headers listed below are also taken to be permitted.
%
% INPUTS:
%   -> none
%
% OUTPUTS:
%   -> PermittedHeaders: structure with:
%       - mnemonics: cell string array of permitted headers
%       - values: cell string array of permitted headers
%       - lowerBounds: cell string array of permitted headers
%       - upperBounds: cell string array of permitted headers

%% DEFINE CONFIGURATION
PermittedHeaders.mnemonics = {'mnemonics','mnems'};
PermittedHeaders.values = {'values','vals'};
PermittedHeaders.lowerBounds = {'lower bounds','lbs'};
PermittedHeaders.upperBounds = {'upper bounds','ubs'};

end

%% FUNCTION TO SETUP FALSE LOGICALS FOR EACH PARAMETER TYPE 
function HeaderAlreadyFound = setup_header_existence_structure(...
    permittedParamVecTypes)
% This function initialises a header existence structure to false.
%
% INPUTS:
%   -> permittedParamVecTypes: cell string array of vector types
%
% OUTPUTS:
%   -> HeaderAlreadyFound: structure with:
%       - values: false
%       - lowerBounds: false
%       - upperBounds: false

%% SETUP STRUCTURE
nPermittedTypes = size(permittedParamVecTypes,2);
for iType = 1:nPermittedTypes
    HeaderAlreadyFound.(permittedParamVecTypes{iType}) = false;
end

end