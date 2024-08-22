function modelFileConfig = get_NLBL_model_file_config 
% model_file_configuration contains information about the format of .maps. 
% It contains information about the layout of *.maps model files. In
% particular, it details the file keyword for each model component. In 
% addition, it details the file layout of information under each of the
% components. It also details whether the particular model component is
% compulsory and which of the components in the file layout are compulsory.
%
% INPUTS:   
%   -> None
%
% OUTPUTS:  
%   -> modelFileConfig - a cell with the following columns:
%       - model file keyword for the object (eg. METADATA)
%       - compulsory flag for that object
%       - structure of & MAPS names for the info under that object
%       - compulsory flags for each of those bits of info
%
% CALLS:
%   -> None
%
% DETAILS:
%
% NOTES:
%
% This version: 29/01/2013
% Author(s): Matt Waldron & Kate Reinold

%% LIST MODEL FILE CONFIGURATION
modelFileConfig = {...
    'METADATA'                     'compulsory'     'metadataFields:metadataDescriptors'                'compulsory:compulsory';
    'ENDOGENOUS VARIABLES'         'compulsory'     'yNames:yMnems'                                     'compulsory:compulsory';
    'EXOGENOUS VARIABLES'          'optional'       'xNames:xMnems'                                     'compulsory:compulsory';
    'RESIDUALS'                    'optional'       'zNames:zMnems'                                     'compulsory:compulsory';
    'PARAMETERS'                   'optional'       'thetaNames:thetaMnems:theta'                       'compulsory:compulsory:optional';
    'EQUATIONS'                    'compulsory'     'yEqNames:yEqMnems'                                 'compulsory:compulsory'};