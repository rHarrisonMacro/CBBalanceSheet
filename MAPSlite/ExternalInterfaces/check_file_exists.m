function check_file_exists(fileName,fileExtension)
% This helper function checks for the existence of a particular file.
% It throws an error if the file specified on input does not exist or if it
% is not of the specified type.
%
% INPUTS:   
%   -> fileName: string name of the file
%   -> fileExtension (optional): type of file - eg 'mat' or 'xlsx'
%
% OUTPUTS:  
%   -> none
%
% DETAILS:  
%   -> This helper validates the existence of any file.
%   -> The input file name can be a full path file name, partial path file
%      name or just the name of a file sitting on the current directory
%      consistent with the logic in the MATLAB "exist" function or in any 
%      load or save commands.
%   -> The check can be optionally file type specific. If the check is file
%      type specific, then the 2nd "fileExtension" input must be passed in. 
%      There are two ways in which the file type is checked: 
%           a) the input file name is specified with an extension, but no
%              specific file extension is input.
%           b) the input file name is specified with an extension and a
%              specific file extension is input (allowing for 
%              differentiation between files of the same name, but not type
%
% NOTES:
%   -> This helper is used in the validation of any MAPS functions that 
%      perform a loading operation. 
%
% This version: 05/12/2013
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
end

%% CHECK FILE NAME IS VALID
if nargin > 1
    check_file_name_is_valid(fileName,fileExtension)
else
    check_file_name_is_valid(fileName);
end

%% CHECK FILE EXISTS
if exist(fileName,'file') ~= 2
    errId = ['MAPS:',mfilename,':NonExistentFile'];
    printableFileName = create_printable_file_name_string(fileName);
    generate_and_throw_MAPS_exception(errId,{printableFileName});
end
    
end