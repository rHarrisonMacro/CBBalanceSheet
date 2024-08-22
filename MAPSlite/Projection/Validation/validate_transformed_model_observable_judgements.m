function validate_transformed_model_observable_judgements(...
    Yj,YjMetadata,SYjs)
% This helper validates an incomplete, transformed observable dataset.
% If any non-real (complex numbers with imaginary components) or inf (and
% -inf) or NaN values are found, then it throws an exception which details 
% the cause of the validation failure, the offending series' metadata and 
% an explanation of what might have caused the problem.
%
% INPUTS:
%   -> Yj: matrix of model observable judgement
%   -> YjMetadata: cell string array of model observable metadata
%   -> SYjs: horizon of each individual judgement
%
% OUTPUTS:
%   -> none
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> validate_incomplete_model_dataset
%   -> validate_model_dataset_does_not_contain_nan
%
% DETAILS:
%   -> This forecast macro helper validates that an incomplete model 
%      observable dataset of judgements does not contain non-real or inf
%      (and -inf) or NaN values.
%   -> This ia a useful check to make after the transformation of
%      judgements made in observable space because the combination of the
%      data transformation equations and the raw observable data (and data
%      for any time-varying trends) may result in model observable data
%      which exhibit those characteristics even if the original raw
%      observable judgements (and time-varying trend judgements) do not
%      (eg. log(-1), log(0), 0/0 etc).
%   -> If all the data passed in is real and non-inf (and the number/shape
%      of the inputs is as expected), then this function returns nothing.
%   -> If not, it throws an exception detailing the causes of the
%      validation failure and the possible underlying causes.
%
% NOTES:
%   -> See <> for a description of data errors exception handling in MAPS.
%
% This version: 08/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of the inputs is as expected. Both inputs
% are compulsory.
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isnumeric(Yj) || ndims(Yj)~=2
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~iscellstr(YjMetadata) || size(YjMetadata,2)~=1
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
end

%% CONSTRUCT A MASTER EXCEPTION
% Construct a master exception with a message detailing the possible
% underlying causes of the validation failure to pass into the validation
% function.
errId = ['MAPS:',mfilename,':BadObsData'];
BadObsJudgementE = generate_MAPS_exception(errId);

%% CHECK DATASET IS REAL & DOES NOT CONTAIN INFS
% Call a routine to check that the dataset does not contain complex numbers
% with imaginary components or inf (or -inf valies). This function will 
% return an exception regardless of the outcome of the validation.
BadObsJudgementE = validate_incomplete_model_dataset(...
    Yj,YjMetadata,BadObsJudgementE);

%% CHECK DATASET DOES NOT CONTAIN NANS
% Extend the incomplete dataset routines for validation that the dataset
% does not contain nans. Catch any exceptions found and add to the master
% exception. Note that this validation allows for the fact that each
% individual judgement may contain a continuous set of NaNs at the end of
% the judgement horizon reflecting non-judgement imposition in those
% periods (whereas discontinuous judgements with NaNs as ruled out in
% earier validation).
nYj = size(Yj,1);
for iYj = 1:nYj    
    try
        validate_model_dataset_does_not_contain_nan(...
            Yj(iYj,1:SYjs(iYj)),YjMetadata(iYj));
    catch NanE
        BadObsJudgementE = addCause(BadObsJudgementE,NanE);
    end
end

%% THROW ANY EXCEPTIONS CAUGHT
% If the master exception contains any causes, throw it back out of this
% function.
if ~isempty(BadObsJudgementE.cause)
    throw(BadObsJudgementE);
end

end