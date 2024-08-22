function validate_retransformed_raw_observable_dataset(Model,Ytilde)
% This helper validates a complete, retransformed raw observable dataset.
% If any non-real (complex numbers with imaginary components), inf (and 
% -inf) or nan values are found, then it throws an exception which details
% the cause of the validation failure, the offending series' metadata and
% an explanation of what might have caused the problem.
%
% INPUTS:   
%   -> Model: MAPS model structure 
%   -> Ytilde: matrix of raw observable data
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> validate_complete_model_dataset
%
% DETAILS:  
%   -> This forecast macro helper validates that a complete raw observable 
%      dataset does not contain non-real, inf (and -inf) or nan values. 
%   -> This ia a useful check to make after the transformation of a 
%      complete set of raw observable data to model observable space 
%      because the combination of the data transformation equations and the 
%      raw observable data (and data for any time-varying trends) may 
%      result in model observable data which exhibit those characteristics 
%      even if the original raw observable data (and time-varying trend 
%      data) does not (eg. log(0), log(-1) etc).
%   -> If all the data passed in is real and finite (and the number / shape 
%      of the inputs is as expected), then this function returns no output.
%   -> If not, it throws an exception detailing the causes of the 
%      validation failure and the possible underlying causes.    
%
% NOTES:   
%   -> See <> for a description of data errors exception handling in MAPS.
%
% This version: 22/02/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of the inputs is as expected. Both inputs 
% are compulsory.
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});     
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~isnumeric(Ytilde) || ndims(Ytilde)~=2
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
end

%% CONSTRUCT A MASTER EXCEPTION
% Construct a master exception with a message detailing the possible
% underlying causes of the validation failure to pass into the validation
% function.
errId = ['MAPS:',mfilename,':BadObsData'];
BadObsDataE = generate_MAPS_exception(errId);

%% CALL VALIDATION ROUTINE
% Call the validation routine for validation of a complete dataset (ie no
% nans). Pass in the mnemonic metadata name for raw observables and the
% exception constructed above. The routine will add any exception causes
% encountered as causes to the exception passed in and throw it back out.
validate_complete_model_dataset(Model,Ytilde,'YtildeMnems',BadObsDataE);

end