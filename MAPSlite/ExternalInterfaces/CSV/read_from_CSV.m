function dataCellArr = read_from_CSV(fileName)
% This helper can be used to read information from a generic CSV file.
% That information can be any mixture of numerics and strings. Each
% separate piece of information is represented in the output as an element
% in a cell array.
%
% INPUTS:   
%   -> fileName: csv filename
%
% OUTPUTS:  
%   -> dataCellArr: cell array of data from the CSV file
%
% DETAILS:  
%   -> This function will read an arbitrary mix of characters and numerics 
%      from a CSV file. The numerics are converted from strings to numbers,
%      which means they are represented as doubles in MATLAB.
%   -> Each individual line in the CSV file is represented as a separate 
%      row in the output cell array. Information between the comma 
%      delimiters on any individual line is represented in separate columns 
%      in the output cell array.
%
% NOTES:
%   -> This function assumes that all information in the CSV file is
%      either numeric or string. It should not be used with other data 
%      types because they may not be read in correctly.
%   -> It can be used in custom macros to read in additional inputs
%      required for the operation of macros above and beyone those provided
%      from the latest run or the provisional run in EASE.
%
% This version: 20/01/2013
% Author(s): Matt Waldron

%% CHECK INPUT
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~ischar(fileName)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
end

%% ATTEMPT TO OPEN FILE
% Check file exists and if it does attempt to open it, throwing an error if
% there are any problems.
check_file_exists(fileName,'csv');
[fID,errMessage] = fopen(fileName);
if ~isempty(errMessage)
    errId = ['MAPS:',mfilename,':FileOpenFailure'];
    generate_and_throw_MAPS_exception(errId,{errMessage});
end

%% SCAN THE CONTENTS OF THE FILE
% Scan the contents of the file line-by-line. Represent each of those lines
% as a separate element in the cell array.
rawFileContents = textscan(fID,'%s','delimiter','\n');
rawFileContents = strtrim(rawFileContents{1});

%% CLOSE THE FILE
% Close the file. Throw an exception if, for some unspecified reason, 
% MATLAB is unable to.
status = fclose(fID);
if status ~= 0
    errId = ['MAPS:',mfilename,':FileCloseFailure'];
    generate_and_throw_MAPS_exception(errId);
end

%% SPLIT EACH COMMA SEPARATED PIECE OF INFORMATION
% Split each line from the file by the comma delimiter used in CSVs.
% Extract each of those separate lines and represent them in a single cell
% array. Note that it need not be the case that each line has the same
% number of commas. In such cases, where that is not true, this function
% fills up from the first column of the cell array outwards, leaving outer
% columns unfilled as necessary. Note that an alternative implementation
% would require the same number of commas per line in which case an error
% should be thrown.
fileContentsSplit = regexp(rawFileContents,',','split');
nLinesInFile = size(fileContentsSplit,1);
nColsPerLineInFile = cellfun(@(x) size(x,2),fileContentsSplit);
maxColsInFile = max(nColsPerLineInFile);
strCellArr = cell(nLinesInFile,maxColsInFile);
for iLine = 1:nLinesInFile
    iLineCols = nColsPerLineInFile(iLine);
    strCellArr(iLine,:) = [fileContentsSplit{iLine} ...
        repmat({''},[1 maxColsInFile-iLineCols])];
end

%% CONVERT NUMERICS
% Convert any string information that can be represented numerically to
% numbers, leaving character information as strings. Replace empty strings
% with NaN to be consistent with behaviour in MATLAB Excel loading, where
% empty cells are represented as NaNs.
numCellArr = cellfun(@str2num,strCellArr,'UniformOutput',false);
numLogicals = ~cellfun(@isempty,numCellArr);
emptyLogicals = cellfun(@isempty,strCellArr);
dataCellArr = strCellArr;
dataCellArr(numLogicals) = numCellArr(numLogicals);
dataCellArr(emptyLogicals) = {NaN};

end