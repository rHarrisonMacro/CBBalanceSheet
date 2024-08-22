function [directory,fileName,fileExtension] = ...
    split_file_name_into_consituent_parts(fullPathFileName)
% This helper splits a full path file name into is constituent parts.
%
% INPUTS:   
%   -> fullPathFileName: full path file name string
%
% OUTPUTS:  
%   -> directory: directory string
%   -> fileName: file name string
%   -> fileExtension: file extension string
%
% DETAILS:  
%   -> This helper splits a full path file name string into its consituent
%      parts of a directory string, a file name string and a file extension
%      strings. 
%   -> Both the directory and file extenion strings could be empty on 
%      output and the input file name could still be valid because files
%      without extensions are treated as .m files in MATLAB and files can
%      be saved and loaded from the current directory without a full path. 
%   -> This function does not validate the file name. See 
%      check_file_name_is_valid for validation.
%
% NOTES:
%   -> This helper is used in the validation of any MAPS functions that 
%      perform a loading operation.
%   -> Updated to support non-Windows file systems.
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
end

%% EXTRACT DIRECTORY & REMAINDER
comparatorString = ['[^' filesep filesep ']+$'];
[fullPathFileNameSplitByFinalBackSlash,fileNameCell] = regexp(...
    fullPathFileName,comparatorString,'split','match');
directory = fullPathFileNameSplitByFinalBackSlash{1};
if isempty(fileNameCell)
    fileNameWithExtension = '';
else
    fileNameWithExtension = fileNameCell{1};
end

%% EXTRACT FILE EXTENSION & REMAINDER
comparatorString = [filesep '.[^' filesep '.]+$'];
[fileNameSplitByFinalDot,fileExtensionCell] = regexp(...
    fileNameWithExtension,comparatorString,'split','match');
fileName = fileNameSplitByFinalDot{1};
if ~isempty(fileExtensionCell)
    fileExtension = fileExtensionCell{1};
else
    fileExtension = '';
end

end