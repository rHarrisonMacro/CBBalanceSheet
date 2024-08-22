function varargout = load_variables_from_mat_file(matFileName,varargin)
% This helper function loads specified variables from a .mat file.
%
% INPUTS:   
%   -> matFileName: string name of the .mat file
%   -> varargin: list of variable names
%
% OUTPUTS:  
%   -> varargout: variables from the mat file
%
% DETAILS:  
%   -> This helper outputs specified variables from any '.mat' file.
%   -> It will throw an error if the matFileName specified does not exist
%      (or is not of the correct type (i.e. with a '.mat' extension)) or if
%      any of the specified variables do not form part of mat file.
%   -> In particular, any variables specified must have been saved in the
%      mat file as separate variables. For instance, if you had a structure
%      with fields 'EmmaKate' and 'AbiHaddow' then the save command would 
%      have had to have used the '-struct' keyword for those two fields to
%      be extracted as separate variables in this function.
%
% NOTES:
%   -> This helper is used to load particular variables from a mat file. 
%
% This version: 04/10/2012
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~ischar(matFileName)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~is_row_cell_string_array(varargin)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);  
elseif nargin-nargout~=1
    errId = ['MAPS:',mfilename,':BadNarginNargout'];
    generate_and_throw_MAPS_exception(errId);
end

%% LOAD MAT FILE CONTENT AS A STRUCTURE
MatFileContent = load_mat_file_as_structure(matFileName);

%% EXTRACT INDIVIDUAL VARIABLES
nVarsToExtract = size(varargin,2);
isVarMissing = true(nVarsToExtract,1);
varargout = cell(1,nVarsToExtract);
for iVar = 1:nVarsToExtract
   if isfield(MatFileContent,varargin{iVar})
       varargout{iVar} = MatFileContent.(varargin{iVar});
       isVarMissing(iVar) = false;
   end
end

%% THROW EXCEPTION FOR MISSING VARIABLES
if any(isVarMissing)
    errId = ['MAPS:',mfilename,':MissingVariables'];
    generate_MAPS_exception_add_causes_and_throw(...
        errId,varargin',isVarMissing);
end

end