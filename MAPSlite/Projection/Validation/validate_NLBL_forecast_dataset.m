function validate_NLBL_forecast_dataset(Model,EvaluatedData,...
    NLBLdatasetValidationE)
% This helper validates a complete NLBL model forecast run dataset.
% If any non-real (complex numbers with imaginary components), inf (and 
% -inf) or nan values are found, then it throws an exception which details
% the cause of the validation failure, the offending series' metadata and
% an explanation of what might have caused the problem.
%
% INPUTS:   
%   -> Model: MAPS model structure 
%   -> EvaluatedData: Structure containing evaluated data.
%   -> NLBLdatasetValidationE: Master exception providing details of the
%      context in which any dataset validation errors have occured (ie
%      projection, inversion over past etc.).
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> unpack_model
%   -> generate_and_throw_MAPS_exception
%   -> validate_complete_model_dataset
%
% DETAILS:  
%   -> This forecast macro helper validates that a complete NLBL model 
%      dataset does not contain non-real, inf (and -inf) or nan values. 
%   -> If all the data passed in is real and finite (and the number / shape 
%      of the inputs is as expected), then this function returns no output.
%   -> If not, it throws an exception detailing the causes of the 
%      validation failure and the possible underlying causes.    
%
% NOTES:   
%   -> See <> for a description of data errors exception handling in MAPS.
%
% This version: 22/02/2011
% Author(s): Alex Haberis and Matt Waldron 06/06/2011

%% CHECK INPUTS
% Check that the number and shape of the inputs is as expected. Both inputs 
% are compulsory.
if nargin < 3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});     
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~isstruct(EvaluatedData)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~strcmp(class(NLBLdatasetValidationE),'MException')
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
end

%% UNPACK MODEL CHARACTERISTICS
% Call the unpacker to find out which variables types are present in the
% model.
[modelHasEndogenousIdentityVariables,...
    modelHasEndogenousNonIdentityVariables,modelHasResiduals] = ...
     unpack_model(Model,{'modelHasEndogenousIdentityVariables',...
     'modelHasEndogenousNonIdentityVariables','modelHasResiduals'});

%% VALIDATE COMPLETE MODEL DATA SET
% Call the validation routine for validation of a complete residuals 
% dataset. Pass in the mnemonic metadata name for residuals and the
% master exception. The routine will add any exception causes
% encountered as causes to the exception passed in and throw it back out.
if modelHasResiduals
    z = EvaluatedData.Residuals;
    try
        validate_complete_model_dataset(Model,z,'zMnems',...
            NLBLdatasetValidationE)
    catch ValidationE
        NLBLdatasetValidationE = ...
            addCause(NLBLdatasetValidationE,ValidationE);
    end
end

%% VALIDATE IDENTITY VARIABLES DATASET
% Call the validation routine for validation of a complete identity
% variables dataset. Pass in the mnemonic metadata name for identity
% variables and the master exception. The routine will add any exception 
% causes encountered as causes to the exception passed in and throw it back
% out.
if modelHasEndogenousIdentityVariables
    yI = EvaluatedData.Identity;
    try
        validate_complete_model_dataset(Model,yI,'yImnems',...
            NLBLdatasetValidationE)
    catch ValidationE
        NLBLdatasetValidationE = ...
            addCause(NLBLdatasetValidationE,ValidationE);
    end
end

%% VALIDATE ENDOGENOUS NON-IDENTITY VARIABLES DATASET
% Call the validation routine for validation of a complete non-identity
% variables dataset. Pass in the mnemonic metadata name for non-identity
% variables and the master exception. The routine will add any exception 
% causes encountered as causes to the exception passed in and throw it back
% out.
if modelHasEndogenousNonIdentityVariables
     yNI = EvaluatedData.NonIdentity;
    try
        validate_complete_model_dataset(Model,yNI,'yNImnems',...
            NLBLdatasetValidationE)
    catch ValidationE
        NLBLdatasetValidationE = ...
            addCause(NLBLdatasetValidationE,ValidationE);
    end
end

%% THROW EXCEPTION
% Throw the exception out of this function if any causes were added.
if ~isempty(NLBLdatasetValidationE.cause)
    throw(NLBLdatasetValidationE);
end

end