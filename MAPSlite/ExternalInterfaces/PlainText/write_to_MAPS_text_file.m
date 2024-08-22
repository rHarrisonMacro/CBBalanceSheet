function write_to_MAPS_text_file(fileContents,fileName)
% This helper function writes out text information to MAPS text files.
% It writes the information out line-by-line, inverting the read operation
% in the read from MAPS text file function.
%
% INPUTS:   
%   -> fileContents: column cell string array
%   -> fileName: full path string name of a *.maps text file
%
% OUTPUTS:  
%   -> none
%
% DETAILS:  
%   -> This helper function writes out string (char) information
%      in the input cell string array to a MAPS text file.
%   -> Each seperate element in the text file is written out to a separate
%      line in the file.
%   -> The operation in this function is consistent with the MAPS text file
%      reader function, which reads in information line-by-line from a MAPS
%      text file.
%
% NOTES:
%   -> One of the applications of this helper is in the creation of MAPS 
%      model files.
%
% This version: 21/02/2013
% Author(s): Matt Waldron

%% CHECK INPUT
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~is_column_cell_string_array(fileContents)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~ischar(fileName)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
end

%% CHECK VALIDITY OF FILE NAME (INC .MAPS EXTENSION)
try
    check_file_name_is_valid(fileName,'maps');
catch InvalidFileNameE
    errId = ['MAPS:',mfilename,':InvalidFileName'];
    generate_MAPS_exception_add_cause_and_throw(...
        InvalidFileNameE,errId,{fileName});
end

%% CALL GENERIC TEXT FILE WRITER
write_to_text_file(fileName,fileContents);

end