function write_to_text_file(txtFileName,contentToWrite,overwriteExisting)
% This helper function writes out to a *.txt or *.maps text file.
%
% INPUTS:   
%   -> txtFileName: full path file name of text file
%   -> contentToWrite: string or column cell string array to write out
%   -> overwriteExisting (optional): true/false
%
% OUTPUTS:  
%   -> none (writes out a file)
%
% DETAILS:  
%   -> This external interface function writes out text (in the form of a
%      set of strings) to a specified text file.
%   -> Each seperate string in the set (stored as a column array) is 
%      written out as a separate line in the text file.
%   -> This function can be used to either: a) append information to a
%      pre-existing text file; b) overwrite a pre-existing text file. In
%      both cases if the text file does not already exist, it is created.
%
% NOTES:
%   -> This is part of a library of MAPS external interface functions.
%   -> If the caller wishes to append or overwrite and existing text file,
%      then the text file name passed in must include the extension (either
%      '.txt' or '.maps' which uniquely identifies that file (or if it does
%      not include the extension, then the file to be amended must have the
%      default extension).
%
% This version: 15/02/2013
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~ischar(txtFileName)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~is_string_or_column_cell_string_array(contentToWrite)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif nargin>2 && ~is_logical_scalar(overwriteExisting)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
end

%% SET DEFAULT FOR OPTIONAL INPUT FLAG
if nargin < 3
    overwriteExisting = true;
end

%% CONVERT CONTENT TO WRITE TO COMMON FORMAT
contentToWrite = convert_string_or_vector_string_array_to_column_array(...
    contentToWrite);

%% CHECK FILE NAME IS VALID
% Check file name is valid and then check that the file extension is 
% either '.txt' or '.maps'. If none was provided, then set the file 
% extension to the default, which is '.txt'.
check_file_name_is_valid(txtFileName);
[~,~,fileExtension] = split_file_name_into_consituent_parts(txtFileName);
if ~isempty(fileExtension)
    validFileExtensions = {'.txt','.maps'};
    if ~any(strcmp(fileExtension,validFileExtensions))
        errId = ['MAPS:',mfilename,':BadFileNameExtension'];
        generate_and_throw_MAPS_exception(...
            errId,{fileExtension validFileExtensions});
    end    
else
    txtFileName = [txtFileName,'.txt'];
end

%% OPEN FILE
% Open the file either as overwrite (w+) or as append (a+) depending on the
% flag. Throw any errors encountered.
if overwriteExisting
    [fID,errMessage] = fopen(txtFileName,'w+');
else
    [fID,errMessage] = fopen(txtFileName,'a+');
end
if ~isempty(errMessage)
    errId = ['MAPS:',mfilename,':FileOpenFailure'];
    generate_and_throw_MAPS_exception(errId,{errMessage});
end

%% WRITE OUT CONTENT
% Compute the number of lines in the content cell array and then write it
% out line-by-line. Capture any exceptions encountered and throw the 
% result. Note the use of \r\n which begins a new line after each string is
% printed (with \r required to ensure a new line in some Microsoft text
% files.
nLinesToWrite = size(contentToWrite,1);
for iLine = 1:nLinesToWrite
    try
        fprintf(fID,'%s\r\n',contentToWrite{iLine});
    catch MATLABwriteE 
        errId = ['MAPS:',mfilename,':FileWriteFailure'];
        generate_MAPS_exception_add_cause_and_throw(...
            MATLABwriteE,errId,contentToWrite(iLine));
    end
end

%% CLOSE THE TEXT FILE
% Close the text file. Throw an error if, for some unspecified reason, 
% MATLAB is unable to.
status = fclose(fID);
if status ~= 0
    errId = ['MAPS:',mfilename,':FileCloseFailure'];
    generate_and_throw_MAPS_exception(errId);
end

end