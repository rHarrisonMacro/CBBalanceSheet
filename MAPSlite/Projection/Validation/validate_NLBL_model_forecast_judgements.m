function validate_NLBL_model_forecast_judgements(JudgementsInfo)
% This NLBL model macro helper validates a forecast judgements dataset.
% It validates the content of the judgements dataset, validates it against
% the model input and checks that the judgements specified meet the rules
% of imposing judgement in NLBL models.
%
% INPUTS:   
%   -> JudgementsInfo: structure with info about all judgements to be 
%      imposed
%
% OUTPUTS:  
%   -> None
%
% CALLS:
%   -> generate_MAPS_exception
%   -> check_NLBL_model_judgement_rule_1 (sub-function)
%   -> check_NLBL_model_judgement_rule_2 (sub-function)
%
% DETAILS:  
%   -> This macro checks NLBL model forecast judgements against their
%      rules:
%        - Rule 1 - the number of residual usages must be equal to the
%          number of endogenous variables that have been fixed.
%        - Rule 2 - it is not possible to fix a residual directly and use
%          it to fix an endogenous variable at the same time.
%   -> If any judgements are found to break the rules, then the macro
%      throws an exception.
%         
% NOTES:
%   -> See <> for a description of MAPS forecast macro helpers.
%   -> See also <> for a list and discussion of the impose judgement rules.
%
% This version: 06/06/2011
% Author(s): Alex Haberis and Matt Waldron 

%% CHECK INPUT
% Complete a basic check on the number and type of inputs passed in.
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)})
elseif ~isstruct(JudgementsInfo)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId)
end

%% CREATE MASTER EXCEPTION TO ADD EACH EXCEPTION ENCOUNTERED
% Create a master exception for the impose judgement rules to which
% exceptions can be added as they are encountered.
masterErrId = ['MAPS:',mfilename,':RulesValidationFailure'];
ImposeJudgementRulesE = generate_MAPS_exception(masterErrId);

%% CHECK NLBL JUDGEMENT RULE 1: NO. RESIDUAL USAGES EQUALS NO. VARS FIXED
% Call sub-function to check NLBL model forecast judgement rule 1.
try 
    check_NLBL_model_judgement_rule_1(JudgementsInfo)
catch Rule1E
    ImposeJudgementRulesE = addCause(ImposeJudgementRulesE,Rule1E);
end

%% CHECK NLBL JUDGEMENT RULE 2: CANNOT FIX RESIDUAL AND USE FOR INVERSION
% Call sub-function to check NLBL model forecast judgement rule 2.
try 
    check_NLBL_model_judgement_rule_2(JudgementsInfo)
catch Rule2E
    ImposeJudgementRulesE = addCause(ImposeJudgementRulesE,Rule2E);
end

%% THROW THE MASTER EXCEPTION IF ANY RULES WERE BROKEN
if ~isempty(ImposeJudgementRulesE.cause)
    throw(ImposeJudgementRulesE)
end
    
end


%% CHECK NLBL JUDGEMENT RULE 1: NO. RESIDUAL USAGES EQUALS NO. VARS FIXED
function check_NLBL_model_judgement_rule_1(JudgementsInfo)
% This NLBL model macro helper validates a forecast judgements dataset.
% It validates the content of the judgements dataset against NLBL model 
% forecast rule 1: the number of residual usages must be equal to the
% number of endogenous variables that have been fixed.
%
% INPUTS:   
%   -> JudgementsInfo: structure with info about all judgements to be 
%      imposed
%
% OUTPUTS:  
%   -> None
%
% CALLS:
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause
%   -> convert_numeric_column_vector_to_string_equivalent
%   -> create_comma_separated_list

%% GENERATE MASTER EXCEPTION
masterErrId = ['MAPS:',mfilename,':Rule1ValidationFailure'];
Rule1E = generate_MAPS_exception(masterErrId);

%% CARRY OUT RULE CHECK
% Check in the JudgementsInfo structure that the number of residual usages
% is neither greater than or less than the the number of variables fixed.
% Since the combination of fixes may change in each of the H quarters of
% the forecast, the check is carried out period by period.  If there are
% any quarters for which the number of residual usages and variables fixed
% is different, then an exception is generated with details of the invalid
% judgements.
if isfield(JudgementsInfo,'residualUsage')
    residualUsage = JudgementsInfo.residualUsage;
    varsFixed = JudgementsInfo.varsFixed;
    H = JudgementsInfo.H;
    nResidualUsages = zeros(1,H);
    nVarsFixed = zeros(1,H);
    for iPeriod = 1:H
        nResidualUsages(iPeriod) = size(residualUsage{iPeriod},1);
        nVarsFixed(iPeriod) = size(varsFixed{iPeriod},1);
    end
    numResUsagesExceedsFixesInds = (nResidualUsages>nVarsFixed);            %#ok<NASGU>
    numFixesExceedsResUsagesInds = (nVarsFixed>nResidualUsages);            %#ok<NASGU>
    indTypes = {'numResUsagesExceedsFixes';...
        'numFixesExceedsResUsages'};
    nIndTypes = size(indTypes,1);
    for iIndType = 1:nIndTypes
        iInds = eval([indTypes{iIndType},'Inds']);
         if any(iInds)
             errId = [masterErrId,':',indTypes{iIndType}];
             Rule1iE = generate_MAPS_exception(errId);
             errId = [masterErrId,':Instance'];
             forecastQuarters = 1:H;
             rule1InvalidQrtrs = forecastQuarters(iInds);
             rule1InvalidQrtrsStrs = ...
                convert_numeric_column_vector_to_string_equivalent(...
                rule1InvalidQrtrs');
             errArgs = {create_comma_separated_list(rule1InvalidQrtrsStrs)};
             Rule1iE = generate_MAPS_exception_and_add_as_cause(...
                    Rule1iE,errId,errArgs);
            Rule1E = addCause(Rule1E,Rule1iE);
         end
             
    end
end

%% THROW THE MASTER EXCEPTION IF THE RULE WAS BROKEN
% If the master exception contains a cause, throw it out of this function.
if ~isempty(Rule1E.cause)
    throw(Rule1E)
end

end


%% CHECK NLBL JUDGEMENT RULE 2: CANNOT FIX RESIDUAL AND USE FOR INVERSION
function check_NLBL_model_judgement_rule_2(JudgementsInfo)
% This NLBL model macro helper validates a forecast judgements dataset.
% It validates the content of the judgements dataset against NLBL model 
% forecast rule 2: it is not possible to fix a residual directly and use it
% in an inversion (ie for an endogenous variable fix).
%
% INPUTS:   
%   -> JudgementsInfo: structure with info about all judgements to be 
%      imposed
%
% OUTPUTS:  
%   -> None
%
% CALLS:
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause
%   -> convert_numeric_column_vector_to_string_equivalent
%   -> create_comma_separated_list

%% GENERATE MASTER EXCEPTION
masterErrId = ['MAPS:',mfilename,':Rule2ValidationFailure'];
Rule2E = generate_MAPS_exception(masterErrId);

%% CARRY OUT CHECK
% Check that a residual hasn't been fixed and used to fix an endogenous
% variable at the same time.  Since the combination of fixes can change in
% each of the H quarters of the forecast, the check is carried out for each
% period.  If the rule is broken at any point, an excelption is generated
% with details of the invalid judgement.
if isfield(JudgementsInfo,'resFixed') && ...
        isfield(JudgementsInfo,'residualUsage')
    resFixed = JudgementsInfo.resFixed;
    residualUsage = JudgementsInfo.residualUsage;
    H = JudgementsInfo.H;
    rule2InvalidFixesInd = zeros(1,H);
    for iPeriod = 1:H
        if ~isempty(resFixed{iPeriod}) && ~isempty(residualUsage{iPeriod})
            rule2InvalidFixesInd(iPeriod) = ...
                any(ismember(resFixed{iPeriod},residualUsage{iPeriod}));
        else
            rule2InvalidFixesInd(iPeriod) = false;
        end
    end
    rule2InvalidFixes = (rule2InvalidFixesInd==1);
	if any(rule2InvalidFixes)
        errId = [masterErrId,':Instance'];
        forecastQuarters = 1:H;
        rule2InvalidQrtrs = forecastQuarters(rule2InvalidFixes);
        rule2InvalidQrtrsStrs = ...
            convert_numeric_column_vector_to_string_equivalent(...
            rule2InvalidQrtrs');
        errArgs = {create_comma_separated_list(rule2InvalidQrtrsStrs)};
        Rule2Ei = generate_MAPS_exception(errId,errArgs);
        Rule2E = addCause(Rule2E,Rule2Ei);
	end
end

%% THROW THE MASTER EXCEPTION IF THE RULE WAS BROKEN
% If the master exception contains a cause, throw it out of this function.
if ~isempty(Rule2E.cause)
    throw(Rule2E)
end
end
