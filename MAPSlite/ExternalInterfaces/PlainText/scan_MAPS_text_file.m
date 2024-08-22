function [fileContents,fileLineNumbers] = scan_MAPS_text_file(fileName)
% This helper function scans in text info from MAPS text files.
% It scans in the text, trims it and removes empty lines / comments.
%
% INPUTS:   
%   -> fileName: full path string name of the *.maps text file
%
% OUTPUTS:  
%   -> fileContents: column cell array with all lines in the file
%   -> fileLineNumbers: line numbers of the file corresponding to each
%      element of fileContents
%
% CALLS:    
%   -> generate_and_throw_MAPS_exception
%
% DETAILS:  
%   -> scan_MAPS_text_file reads in information from MAPS related text 
%      files. 
%   -> It trims that information for leading and trailing blank spaces, 
%      empty lines and commented out information (as inferred from the % 
%      symbol).
%
% NOTES:
%   -> This helper is used in the parsing of MAPS models. See <> for 
%      information about the format of MAPS model files.
%
% This version: 11/01/2011
% Author(s): Matt Waldron

%% CHECK INPUT
% Check that the number and type of input is as expected by the this 
% function.
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~ischar(fileName)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
end

%% ATTEMPT TO OPEN FILE
% Attempt to open the model file represented by the model file name. Throw
% an error if the file doesn't have the correct *.maps appendage, if it 
% does not exist or if its not possible to open the file.
if isempty(strfind(fileName,'.maps'))
    errId = ['MAPS:',mfilename,':MissingFileNameAppendage'];
    generate_and_throw_MAPS_exception(errId);
elseif exist(fileName,'file') ~= 2
    errId = ['MAPS:',mfilename,':UnknownFile'];
    generate_and_throw_MAPS_exception(errId);
end
[fID,errMessage] = fopen(fileName);
if ~isempty(errMessage)
    errId = ['MAPS:',mfilename,':FileOpenFailure'];
    generate_and_throw_MAPS_exception(errId,{errMessage});
end

%% SCAN IN CONTENT OF TEXT FILE
% Scan in the content of the text file with the text scan delimiter and
% extract the information as a single cell array.
rawFileContents = textscan(fID,'%s','delimiter','\n');
rawFileContents = strtrim(rawFileContents{1});

%% CLOSE THE TEXT FILE
% Close the text file. Throw an error if, for some unspecified reason, 
% MATLAB is unable to.
status = fclose(fID);
if status ~= 0
    errId = ['MAPS:',mfilename,':FileCloseFailure'];
    generate_and_throw_MAPS_exception(errId);
end

%% REMOVE COMMENTED OUT INFORMATION
% Search for the % comment delimiter in the cell array. Remove all 
% information found after the first delimiter in that line.
rawFileContentsSplit = regexp(rawFileContents,'%','split');
nLinesInFile = size(rawFileContents,1);
for iLineInFile = 1:nLinesInFile
    rawFileContents{iLineInFile} = strtrim(...
        rawFileContentsSplit{iLineInFile}{1});
end

%% REMOVE BLANK LINES
% Remove any blank lines. (This will also remove entrie lines that were 
% commented out as split in the cell above). If all the lines were blank or 
% the cell array is of an unexpected shape, throw an error.
nLinesInFile = size(rawFileContents,1);
linesInFile = (1:nLinesInFile)';
emptyLineLogicals = cellfun(@isempty,rawFileContents);
fileContents = rawFileContents(~emptyLineLogicals);
fileLineNumbers = linesInFile(~emptyLineLogicals);
if isempty(rawFileContents)
    errId = ['MAPS:',mfilename,':EmptyScannedContent'];
    generate_and_throw_MAPS_exception(errId);
end
if (size(rawFileContents,1)<2) && (size(rawFileContents,2)~=1)
    errId = ['MAPS:',mfilename,':UnexpectedScannedContent'];
    generate_and_throw_MAPS_exception(errId);
end

end