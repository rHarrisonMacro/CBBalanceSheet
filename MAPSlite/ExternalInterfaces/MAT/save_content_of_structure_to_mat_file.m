function save_content_of_structure_to_mat_file(matFileName,StructToSave)
% This helper function saves the content (fields) of a structure to .mat.
%
% INPUTS:   
%   -> matFileName: string name of the .mat file
%   -> StructToSave: structure whose contents to save out
%
% OUTPUTS:  
%   -> none (saves a .mat file)
%
% DETAILS:  
%   -> This helper saves the content of a structure to .mat.
%   -> It invokes the '-struct' option in the save command, which means 
%      that the fields of the structure are saved out as separate
%      variables.
%   -> The structure can be loaded back in using another MAPS .mat
%      interface helper.
%   -> It will throw an error if the matFileName specified is not valid.
%      Note that it does not require the filename to have a .mat extension
%      because that is the default file type in MATLAB.
%
% NOTES:
%   -> This helper is part of a family of .mat interface functions in MAPS. 
%
% This version: 22/02/2013
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~ischar(matFileName)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(StructToSave)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);       
end

%% CHECK VALIDITY OF FILE NAME
try
    check_file_name_is_valid(matFileName);
catch FileNameE
    errId = ['MAPS:',mfilename,':BadFileName'];
    generate_MAPS_exception_add_cause_and_throw(...
        FileNameE,errId,{matFileName});
end

%% ATTEMPT TO SAVE THE CONTENT OF THE STRUCTURE
% The content of the structure is saved with each field of the structure as
% a separate variable in the resulting .mat file - this is enforced by the 
% "-struct" option in the call to save.
try    
    save(matFileName,'-struct','StructToSave');
catch SaveE
    errId = ['MAPS:',mfilename,':UnableToSaveStruct'];
    generate_MAPS_exception_add_cause_and_throw(...
        SaveE,errId,{matFileName});
end

end