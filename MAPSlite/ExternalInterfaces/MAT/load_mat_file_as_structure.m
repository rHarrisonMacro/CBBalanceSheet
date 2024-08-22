function MatFileContent = load_mat_file_as_structure(matFileName)
% This helper function loads the content of a .mat file as a structure.
%
% INPUTS:   
%   -> matFileName: string name of the .mat file
%
% OUTPUTS:  
%   -> MatFileContent: structure containing the content of the .mat file
%
% CALLS:    
%   -> generate_and_throw_MAPS_exception
%   -> check_file_exists
%
% DETAILS:  
%   -> This helper loads the content of any '.mat' file as a structure.
%   -> It will throw an error if the matFileName specified does not exist
%      or is not of the correct type (i.e. with a '.mat' extension.
%
% NOTES:
%   -> This helper is used to check and load '.mat' files. 
%
% This version: 04/10/2012
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~ischar(matFileName)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
end

%% CHECK FOR EXISTENCE OF FILE
check_file_exists(matFileName,'.mat');

%% LOAD MAT FILE CONTENT AS A STRUCTURE
MatFileContent = load(matFileName,'-mat');

end