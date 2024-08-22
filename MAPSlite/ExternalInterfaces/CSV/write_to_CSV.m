function write_to_CSV(dataCellArr,fileName)
% This helper can be used to write generic information to a CSV file.
% That information can be any mixture of numerics and strings. Each
% separate piece of information must be represented in the data input as an 
% element in a cell array.
%
% INPUTS:
%   -> dataCellArr: cell array of data
%   -> fileName: CSV filename
%
% OUTPUTS:  
%   -> none
%
% DETAILS:  
%   -> This function will write an arbitrary mix of characters and numerics 
%      to a CSV file. The numerics are converted to strings (with a 
%      precision of 16 significant figures).
%   -> Each separate row in the input cell array is represented as a 
%      separate line in the CSV file. Information in separate columns in 
%      the cell array is separated by comma delmiters in the CSV file. 
%
% NOTES:
%   -> This function requires that all information in the input cell array
%      be either numeric or string.
%   -> It can be used in custom macros to write additional information as
%      output or diagnostic information on the runs executed.
%
% This version: 20/01/2013
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~is_two_dimensional_cell_array(dataCellArr)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~ischar(fileName)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);    
end

%% CONVERT DATA CELL ARRAY NUMERICS TO STRINGS
% If any of the contents of the data cell array input are not numerics or
% characters, throw an exception. Otherwise, convert the data cell array to
% a string equivalent where all numeric information is represented by
% strings and where NaNs have been replaced by ''.
strLogicals = cellfun(@ischar,dataCellArr);
numLogicals = cellfun(@isnumeric,dataCellArr);
nanLogicals = cellfun(...
    @(x) isnumeric(x)&&isnan(x),dataCellArr,'UniformOutput',true);
if ~all(all(strLogicals|numLogicals))
    errId = ['MAPS:',mfilename,':UnhandledDataType'];
    generate_and_throw_MAPS_exception(errId);
end
strCellArr = dataCellArr;
strCellArr(numLogicals) = cellfun(...
    @(x) num2str(x,16),strCellArr(numLogicals),'UniformOutput',false);
strCellArr(nanLogicals) = {''};

%% ADD COMMA DELIMITERS
% Add comma delimiters to the cell array to sepaerate each column and then
% concatenate the columns into a one column cell array in preparation for
% writing.
[nLinesInFile,nColsInFile] = size(strCellArr);
rawFileContentsSplit = strcat(strCellArr',...
    [repmat({','},[nColsInFile-1 nLinesInFile]);...
    repmat({''},[1 nLinesInFile])])';
rawFileContents = cell(nLinesInFile,1);
for iCol = 1:nColsInFile
   rawFileContents = strcat(rawFileContents,rawFileContentsSplit(:,iCol));   
end

%% ATTEMPT TO OPEN FILE
% Check file name is valid and if it is attempt to open the file 
% corresponding the file name input. If it could not be opened for some 
% reason then throw an exception.
check_file_name_is_valid(fileName,'csv');
[fID,errMessage] = fopen(fileName,'w+');
if ~isempty(errMessage)
    errId = ['MAPS:',mfilename,':FileOpenFailure'];
    generate_and_throw_MAPS_exception(errId,{errMessage});
end

%% WRITE OUT THE DATA
% Write out the data line-by-line. If for any reason, the print command
% fails, throw an exception detailing the cause of the failure.
for iLine = 1:nLinesInFile
    try
        fprintf(fID,'%s\n',rawFileContents{iLine});
    catch MATLABwriteE 
        errId = ['MAPS:',mfilename,':FileWriteFailure'];
        generate_MAPS_exception_add_cause_and_throw(MATLABwriteE,errId,...
            rawFileContents(iFileLine));
    end
end

%% CLOSE THE FILE
% Close the file. Throw an exception if, for some unspecified reason, 
% MATLAB is unable to.
status = fclose(fID);
if status ~= 0
    errId = ['MAPS:',mfilename,':FileCloseFailure'];
    generate_and_throw_MAPS_exception(errId);
end

end