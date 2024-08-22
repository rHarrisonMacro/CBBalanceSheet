function check_file_name_is_valid(fullPathFileName,expectedExtension)
% This helper checks that a file name is valid.
%
% INPUTS:   
%   -> fullPathFileName: full path file name string
%   -> expectedExtension (optional): type of file - eg 'mat' or 'xlsx'
%
% OUTPUTS:  
%   -> none
%
% DETAILS:  
%   -> This helper checks that a file name is valid. 
%   -> A file name is valid if the directory embedded in it is valid, the
%      filename itself is valid (contains no invalid characters) and if the
%      file extension is valid.
%
% NOTES:
%   -> This helper is used in the validation of any MAPS functions that 
%      perform a loading operation.
%   -> The code has been updated to make it more accessible to alternative
%      filesystems (Windows/Linux/iOS) by using the filesep representation
%      for the file seperator.
%
% This version: 14/12/2022
% Author(s): Matt Waldron & Richard Harrison

%% CHECK INPUTS
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~ischar(fullPathFileName)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif nargin>1 && ~ischar(expectedExtension)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
end

%% CHECK OPTIONAL EXPECTED FILE EXTENSION INPUT
if nargin > 1
    checkFileExtension = true;
    fileSeparationString = [filesep '.'];
    if isempty(regexp(expectedExtension,fileSeparationString,'match'))
        expectedExtension = ['.',expectedExtension];
    end
    errId = ['MAPS:',mfilename,':BadlySpecifiedExpectedExtension'];
    check_file_extension_is_valid(expectedExtension,errId);
else
    checkFileExtension = false;
end

%% SPLIT FILE NAME INTO CONSTITUENT PARTS
[directory,fileName,fileExtension] = ...
    split_file_name_into_consituent_parts(fullPathFileName);

%% CHECK DIRECTORY
if ~isempty(directory) && ~isdir(directory)
    errId = ['MAPS:',mfilename,':BadDirectory'];
    printableDir = create_printable_file_name_string(directory);
    generate_and_throw_MAPS_exception(errId,{printableDir});
end

%% CHECK FILE NAME
[invalidCharsExpr,illegalChars] = ...
    get_regexp_expression_for_invalid_file_name_chars;
if ~isempty(regexp(fileName,invalidCharsExpr,'match'));
    errId = ['MAPS:',mfilename,':BadFileName'];
    printableFileName = create_printable_file_name_string(fileName);
    generate_and_throw_MAPS_exception(...
        errId,{printableFileName illegalChars});
end

%% CHECK FILE EXTENSION
if ~isempty(fileExtension)
    errId = ['MAPS:',mfilename,':BadFileExtension'];
    check_file_extension_is_valid(fileExtension,errId);
end
if checkFileExtension
    if ~strcmp(fileExtension,expectedExtension)
        errId = ['MAPS:',mfilename,':FileExtensionMismatch'];
        generate_and_throw_MAPS_exception(...
            errId,{fileExtension expectedExtension});
    end
end

end

%% HELPER FUNCTION TO CHECK FILE EXTENSION IS VALID
function check_file_extension_is_valid(fileExtension,errId)
% This function compiles the regexp expression for invalid file name chars.
%
% INPUTS:   
%   -> fileExtension: file extension string like ".mat"
%   -> errId: identifier to construct any exception with
%
% OUTPUTS:  
%   -> none
%
% CALLS:    
%   -> none

%% EXTRACT VALID FILE EXTENSION
extensionCheckStr = [filesep '.[a-z|A-Z]+$'];
validFileExtension = regexp(fileExtension,extensionCheckStr,'match');

%% CHECK
if isempty(validFileExtension) || ~strcmp(fileExtension,validFileExtension)
    generate_and_throw_MAPS_exception(errId,{fileExtension});
end

end

%% HELPER FUNCTION FOR A REGEXP EXPRESSION WITH INVALID FILE NAME CHARS
function [invalidCharsExpr,invalidFileNameChars] = ...
    get_regexp_expression_for_invalid_file_name_chars
% This function compiles the regexp expression for invalid file name chars.
%
% INPUTS:   
%   -> none
%
% OUTPUTS:  
%   -> invalidCharsExpr: regexp expression to serach for illegal characters
%   -> invalidFileNameChars: cell array of illegal characters
%
% CALLS:    
%   -> none

%% DEFINE INVALID FILE NAME CHARACTERS
invalidFileNameChars = {'\' '/' '.'};

%% COMPILE REGEXP EXPRESSION
nInvalidChars = size(invalidFileNameChars,2);
% cellOfInvalidCharsAndDelims = [repmat({'\'},[1 nInvalidChars]);...
%     invalidFileNameChars;...
%     [repmat({'|'},[1 nInvalidChars-1]) {''}]];
cellOfInvalidCharsAndDelims = [repmat({filesep},[1 nInvalidChars]);...
    invalidFileNameChars;...
    [repmat({'|'},[1 nInvalidChars-1]) {''}]];
invalidCharsExpr = [cellOfInvalidCharsAndDelims{:}];

end