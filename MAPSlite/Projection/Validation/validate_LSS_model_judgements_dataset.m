function [JudgementsInclude,JudgementsInfo] = ...
    validate_LSS_model_judgements_dataset(Model,Judgements)
% This LSS model macro helper validates a forecast judgements dataset.
% It validates the content of the judgements dataset, validates it against
% the model input and checks that the judgements specified meet the rules
% of imposing judgement in LSS models.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> Judgements: structure with info about all judgements to be imposed
%       - AnticipatedFixes (optional): structure of anticipated fixes
%           - modelVariables: a cell array of metadata & fix values
%           - modelObservables: a cell array of metadata & fix values
%           - rawObservables: a cell array of metadata & fix values
%           - shockUsages: a cell array of metadata & usage indicators
%       - UnanticipatedFixes (optional): structure of unanticipated fixes
%           - modelVariables: a cell array of metadata & fix values
%           - modelObservables: a cell array of metadata & fix values
%           - rawObservables: a cell array of metadata & fix values
%           - shockUsages: a cell array of metadata & usage indicatorse
%       - Shocks (optional): structure detailing judgements over shocks:
%           - anticipated: a cell array of metadata & judgements
%           - unanticipated: a cell array of metadata & judgements
%       - timeVaryingTrends (optional): cell array of metadata & judgements
%
% OUTPUTS:  
%   -> JudgementsInclude: structure describing the judgements content
%       - anticipatedModelVariableFixes: true or false
%       - anticipatedModelObservableFixes: true or false
%       - anticipatedRawObservableFixes: true or false
%       - shocksForAnticipatedFixes: true or false
%       - unanticipatedModelVariableFixes: true or false
%       - unanticipatedModelObservableFixes: true or false
%       - unanticipatedRawObservableFixes: true or false
%       - shocksForUnanticipatedFixes: true or false
%       - anticipatedShockJudgements: true or false
%       - unanticipatedShockJudgements true or false
%       - timeVaryingTrendJudgements: true or false
%   -> JudgementsInfo: new structure with as many fields as indicated true
%      in the judgement content descriptor input, each of which contain:
%       - metadata: judgement metadata (from judgements dataset)
%       - cell: judgement data as a cell array (from judgements dataset)
%       - modelMetadata: complete model metadata for that variable type
%       - indices: model index numbers equivalent to the metadata
%       - matrix: matrix equivalent to the cell vectors
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> unpack_model
%   -> interrogate_judgements_structure (sub-function)
%   -> create_judgements_info_structure (sub-function)
%   -> validate_judgements_against_rules (sub-function)
%
% DETAILS:  
%   -> This helper function valdiates a MAPS LSS model judgement dataset.
%   -> It first checks whether or not the dataset is of the right shape
%      (i.e. that it is a structure with the expected fields and that the
%      cell arrays in the fields have been constructed correctly).
%   -> Next, it creates two structures: one describing the content of the
%      judgement input and another providing more information about the
%      judgement (eg index numbers of the variables over which judgement is
%      being applied in the model. These two structures can be passed as
%      output to this function in the call.
%   -> Finally, it validates the judgements against the rules of imposing
%      judgment in MAPS LSS models. If the judgements fial validation then 
%      this function will throw an exception detailing the cause of the 
%      validation failure. 
%
% NOTES:
%   -> See <> for a description of MAPS forecast macro helpers.
%   -> See also <> for a list and discussion of the impose judgement rules.
%
% This version: 03/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number & shape of inputs is as expected. Both inputs are
% compulsory.
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(Judgements)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
end

%% UNPACK MODEL CLASS INFO
% Unpack model class information (i.e. whether or not the model is linear
% state space). If it is not linear state space, throw an exception.
modelIsLinearStateSpace = unpack_model(Model,{'modelIsLinearStateSpace'});
if ~modelIsLinearStateSpace
    errId = ['MAPS:',mfilename,':BadModelClass'];
    generate_and_throw_MAPS_exception(errId);
end

%% INTERROGATE JUDGEMENT STRUCTURE
% Interrogate the judgement structure to determine its content. This
% function will throw an exception if the content of the judgement
% structure is not as expected or if it is not consistent with the model
% being used.
JudgementsInclude = interrogate_judgements_structure(Model,Judgements);

%% CREATE JUDGEMENT INFO STRUCTURE
% Create a judgement info structure in preparation for validating the
% judgements against the LSS model impose judgement rules. This function
% will throw an exception if it is unable to convert the metadata
% describing which variables to impose judgement on to equivalent model
% index numbers.
JudgementsInfo = create_judgements_info_structure(...
    Model,Judgements,JudgementsInclude);

%% VALIDATE THAT THE JUDGEMENTS MEET THE IMPOSE JUDGEMENT RULES
% Call another sub-function to validate that the judgements meet the rules
% of imposing judgement in MAPS LSS models. If the judgements fail
% validation, this function will throw an exception describing which rule
% was violated and why.
validate_judgements_against_rules(Model,JudgementsInclude,JudgementsInfo);

end

%% FUNCTION TO VALIDATE JUDGEMENTS AGAINST THE IMPOSE JUDGEMENT RULES
function validate_judgements_against_rules(...
    Model,JudgementsInclude,JudgementsInfo)
% This helper manages the checking of judgements against the rules.
% It checks that the combined set of judgements meet the rules of
% imposing judgement in MAPS LSS models. If any of them do not, it throws 
% an exception detailing which rule or rules were broken and why. 
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> JudgementsInclude: structure describing the judgements content
%       - anticipatedModelVariableFixes: true or false
%       - anticipatedModelObservableFixes: true or false
%       - anticipatedRawObservableFixes: true or false
%       - shocksForAnticipatedFixes: true or false
%       - unanticipatedModelVariableFixes: true or false
%       - unanticipatedModelObservableFixes: true or false
%       - unanticipatedRawObservableFixes: true or false
%       - shocksForUnanticipatedFixes: true or false
%       - anticipatedShockJudgements: true or false
%       - unanticipatedShockJudgements true or false
%       - timeVaryingTrendJudgements: true or false
%   -> JudgementsInfo: new structure with as many fields as indicated true
%      in the judgement content descriptor input, each of which contain:
%       - metadata: judgement metadata (from judgements dataset)
%       - cell: judgement data as a cell array (from judgements dataset)
%       - modelMetadata: complete model metadata for that variable type
%       - indices: model index numbers equivalent to the metadata
%       - matrix: matrix equivalent to the cell vectors
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> check_judgement_rule_1 (sub-function)
%   -> check_judgement_rule_2 (sub-function)
%   -> check_judgement_rule_3 (sub-function)
%   -> check_judgement_rule_4 (sub-function)
%   -> check_judgement_rule_5 (sub-function)

%% CREATE MASTER EXCEPTION TO ADD EACH EXCEPTION ENCOUNTERED
% Create a master exception for the impose judgement rules to which
% exceptions can be added as they are encountered.
errId = ['MAPS:',mfilename,':RulesValidationFailure'];
ImposeJudgementRulesE = generate_MAPS_exception(errId);

%% DEFINE THE FUNCTIONS USED TO CHECK THE RULES
% Use function handles to define calls to a set of sub-functions which wil
% validate each impose judgement rule. (Each of these functions takes
% the model, judgemnt content and judgement info structure as inputs.
functionsToCheckRules = {
    @(Model,JudgementsInclude,JudgementsInfo) check_judgement_rule_1(...
        Model,JudgementsInclude,JudgementsInfo);
    @(Model,JudgementsInclude,JudgementsInfo) check_judgement_rule_2(...
        Model,JudgementsInclude,JudgementsInfo);
    @(Model,JudgementsInclude,JudgementsInfo) check_judgement_rule_3(...
        Model,JudgementsInclude,JudgementsInfo);
    @(Model,JudgementsInclude,JudgementsInfo) check_judgement_rule_4(...
        Model,JudgementsInclude,JudgementsInfo);
    @(Model,JudgementsInclude,JudgementsInfo) check_judgement_rule_5(...
        Model,JudgementsInclude,JudgementsInfo)
        };

%% CALL EACH OF THE RULE CHECKING FUNCTIONS
% Run through each of the functions to check. If any of them throw
% exceptions, ass the exceptions as causes to the master exception. 
nFunctionsToCall = size(functionsToCheckRules,1);
for iCall = 1:nFunctionsToCall
   try 
       ifunctionToCall = functionsToCheckRules{iCall};
       ifunctionToCall(Model,JudgementsInclude,JudgementsInfo);
   catch RuleE
       ImposeJudgementRulesE = addCause(ImposeJudgementRulesE,RuleE);
   end
end

%% THROW IMPOSE JUDGEMENT RULES EXCEPTION
% If the master exception created above is non-empty, throw it and this
% function will be terminated with an appropriate message.
if ~isempty(ImposeJudgementRulesE.cause)
    throw(ImposeJudgementRulesE)
end

end

%% FUNCTION TO CHECK JUDGEMENTS AGAINST IMPOSE JUDGEMENT RULE 1
function check_judgement_rule_1(~,JudgementsInclude,JudgementsInfo)
% This helper validates judgements againts rule 1.
% Rule 1 is that the number of each type of fix (anticipated and
% unanticipated) cannot exceed the number of each type of shock used to
% implement that fix in any period (except where there are no fixes, in
% which case there must be no shocks).
%
% INPUTS:  
%   -> JudgementsInclude: structure of indicators for each judgement type
%   -> JudgementsInfo: structure of information about judgements being made
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause

%% CREATE A MASTER EXCEPTION FOR RULE 1
% Create a master exception for rule 1 to add exceptions to as encountered
% below.
masterErrId = ['MAPS:',mfilename,':Rule1ValidationFailure'];
Rule1E = generate_MAPS_exception(masterErrId);

%% CHECK EACH TYPE OF FIX AGAINST THE RULES
% Run through each fix type (anticipated and unaticipated) in turn,
% combining all fixes of that type into one matrix.  Compute the number of
% the fix type in each period and compare that to the number of shock types
% for the inversion. Add exceptions to the master exception as appropriate.
fixTypes = {'Anticipated';'Unanticipated'};
nFixTypes = size(fixTypes,1);
for iType = 1:nFixTypes
    shocksField = ['shocksFor',fixTypes{iType},'Fixes'];
    if JudgementsInclude.(shocksField);
        iMat = JudgementsInfo.(shocksField).matrix;
        fixFields = {[lower(fixTypes{iType}),'ModelVariableFixes'];
            [lower(fixTypes{iType}),'ModelObservableFixes'];
            [lower(fixTypes{iType}),'RawObservableFixes']};
        nFixFields = size(fixFields,1);
        jMat = [];
        for iField = 1:nFixFields
            if JudgementsInclude.(fixFields{iField});
                jMat = [jMat;JudgementsInfo.(fixFields{iField}).matrix];    %#ok<AGROW>
            end
        end
        njPerPeriod = sum(~isnan(jMat),1);
        niPerPeriod = sum(~isnan(iMat),1);
        numFixesExceedsShocksInds = find(njPerPeriod>niPerPeriod);          %#ok<NASGU>
        numShocksExceedsZeroFixesInds = find(...
            (njPerPeriod==0).*(niPerPeriod>njPerPeriod));                   %#ok<NASGU>
        indTypes = {'numFixesExceedsShocks';'numShocksExceedsZeroFixes'};
        nIndTypes = size(indTypes,1);
        for iIndType = 1:nIndTypes
            iInds = eval([indTypes{iIndType},'Inds']);
            if ~isempty(iInds)
                errId = [masterErrId,':',indTypes{iIndType}];
                Rule1iE = generate_MAPS_exception(errId);
                errId = [masterErrId,':Instance'];
                errArgs = {lower(fixTypes{iType}) num2str(iInds)};
                Rule1iE = generate_MAPS_exception_and_add_as_cause(...
                    Rule1iE,errId,errArgs);
                Rule1E = addCause(Rule1E,Rule1iE);
            end
        end
    end
end

%% THROW THE MASTER EXCEPTION IF EITHER RULE WAS BROKEN
% If the master exception contains a cause, throw it out of this function.
if ~isempty(Rule1E.cause)
    throw(Rule1E)
end

end

%% FUNCTION TO CHECK JUDGEMENTS AGAINST IMPOSE JUDGEMENT RULE 2
function check_judgement_rule_2(~,JudgementsInclude,JudgementsInfo)
% This helper validates judgements againts rule 2.
% Rule 2 is that a variable being fixed or shocks being used in a fix under
% a particular assumption (either anticipated or unanticipated) cannot also
% be fixed or used in a fix under the alternative assumption in the same
% period.
%
% INPUTS:
%   -> JudgementsInclude: structure of indicators for each judgement type
%   -> JudgementsInfo: structure of information about judgements being made
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> validate_judgement_datasets_do_not_overlap (sub-function)

%% SETUP MASTER EXCEPTION
% Create a master exception for rule 1 to add exceptions to as encountered
% below.
masterErrId = ['MAPS:',mfilename,':Rule2ValidationFailure'];
Rule2E = generate_MAPS_exception(masterErrId);

%% DEFINE PAIRS OF FIXES TO CHECK
% Define pairs of fixes to check. There are four pairs; three for variables
% and one for shocks.
fixPairsToCheck = {...
    'anticipatedModelVariable'   'unanticipatedModelVariable'
    'anticipatedModelObservable' 'unanticipatedModelObservable'
    'anticipatedRawObservable'   'unanticipatedRawObservable'
    'shocksForAnticipated'       'shocksForUnanticipated'};

%% SEARCH OVER PAIRS OF JUDGEMENTS
% Check each pair one-by-one provided that both members of a pair exist in
% the judgements input. If they do, use sub-function find repetitions to
% return the index numbers (in the model) of repetitions and the periods in
% which each of those repetitions applies. If any repetitions are found,
% call a helper sub-function to create an exception with each instance of
% repetition added as cause.
nPairsToCheck = size(fixPairsToCheck,1);
for iPair = 1:nPairsToCheck
    fixToCheck1 = [fixPairsToCheck{iPair,1},'Fixes'];
    fixToCheck2 = [fixPairsToCheck{iPair,2},'Fixes'];    
    if JudgementsInclude.(fixToCheck1) && JudgementsInclude.(fixToCheck2)
        try
            errId = [masterErrId,':Instance'];
            Rule2Ei = generate_MAPS_exception(...
                errId,{fixToCheck1 fixToCheck2});
            validate_judgement_datasets_do_not_overlap(...
                JudgementsInfo.(fixToCheck1).matrix,...
                JudgementsInfo.(fixToCheck1).metadata,...
                JudgementsInfo.(fixToCheck2).matrix,...
                JudgementsInfo.(fixToCheck2).metadata,Rule2Ei);
        catch RepetitionsE
            Rule2E = addCause(Rule2E,RepetitionsE);
        end
    end
end

%% THROW THE MASTER EXCEPTION IF THE RULE WAS BROKEN
% If the master exception contains a cause, throw it out of this function.
if ~isempty(Rule2E.cause)
    throw(Rule2E)
end

end

%% FUNCTION TO CHECK JUDGEMENTS AGAINST IMPOSE JUDGEMENT RULE 3
function check_judgement_rule_3(Model,JudgementsInclude,JudgementsInfo)
% This helper validates judgements againts rule 3.
% Rule 3 is that a variable being fixed in one space cannot also be fixed 
% in another space in the same period.  
%
% INPUTS:
%   -> Model: MAPS LSS model structure
%   -> JudgementsInclude: structure of indicators for each judgement type
%   -> JudgementsInfo: structure of information about judgements being made
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> transform_indices_from_raw_to_model_observable_space
%   -> transform_indices_from_observable_to_model_var_space
%   -> validate_judgement_datasets_do_not_overlap (sub-function)

%% SETUP A MASTER EXCEPTION
% Create a master exception for rule 1 to add exceptions to as encountered
% below.
masterErrId = ['MAPS:',mfilename,':Rule3ValidationFailure'];
Rule3E = generate_MAPS_exception(masterErrId);

%% DEFINE FIX PAIRS
% Define the variable fix pairs to check. For example, it is not prossible
% to impose judgement on GDP in both raw observable and model observable
% space.
fixPairsToCheck = {
    'anticipatedModelObservable'   'anticipatedRawObservable'
    'unanticipatedModelObservable' 'anticipatedRawObservable'
    'anticipatedModelObservable'   'unanticipatedRawObservable'
    'unanticipatedModelObservable' 'unanticipatedRawObservable'
    'anticipatedModelVariable'     'anticipatedRawObservable'
    'unanticipatedModelVariable'   'anticipatedRawObservable'
    'anticipatedModelVariable'     'unanticipatedRawObservable'
    'unanticipatedModelVariable'   'unanticipatedRawObservable'    
    'anticipatedModelVariable'     'anticipatedModelObservable'
    'unanticipatedModelVariable'   'anticipatedModelObservable'
    'anticipatedModelVariable'     'unanticipatedModelObservable'
    'unanticipatedModelVariable'   'unanticipatedModelObservable'};

%% SEARCH OVER RELEVANT PAIRS OF FIXES
% Run through each of the pairs to check. If both members of the pair exist
% in the judgements input, then search for repetitions in usage in any
% periods. If any are found then add a cause to the exception. If one of
% the comparators in the pairs is a raw observable, then convert the index
% numbers of the raw observable fixes to model observable space. If the 
% other comparator is a model variable then convert the index numbers of 
% the model observable fixes to model variable space. 
nPairsToCheck = size(fixPairsToCheck,1);
for iPair = 1:nPairsToCheck
    fixToCheck1 = [fixPairsToCheck{iPair,1},'Fixes'];
    fixToCheck2 = [fixPairsToCheck{iPair,2},'Fixes'];
    if JudgementsInclude.(fixToCheck1) && JudgementsInclude.(fixToCheck2)
        JudgementsInfoToCheck1 = JudgementsInfo.(fixToCheck1);
        JudgementsInfoToCheck2 = JudgementsInfo.(fixToCheck2);
        if ~isempty(strfind(fixToCheck2,'RawObservable'))
            JudgementsInfoToCheck2.indices = ...
                transform_indices_from_raw_to_model_observable_space(...
                Model,JudgementsInfoToCheck2.indices); 
            JudgementsInfoToCheck2.metadata = ...
                JudgementsInfoToCheck1.modelMetadata(...
                JudgementsInfoToCheck2.indices);
        end
        if ~isempty(strfind(fixToCheck1,'ModelVariable'))
            JudgementsInfoToCheck2.indices = ...
                transform_indices_from_observable_to_model_var_space(...
                Model,JudgementsInfoToCheck2.indices);
            JudgementsInfoToCheck2.metadata = ...
                JudgementsInfoToCheck1.modelMetadata(...
                JudgementsInfoToCheck2.indices);            
        end        
        try
            errId = [masterErrId,':Instance'];
            Rule3Ei = generate_MAPS_exception(...
                errId,{fixToCheck1 fixToCheck2});
            validate_judgement_datasets_do_not_overlap(...
                JudgementsInfoToCheck1.matrix,...
                JudgementsInfoToCheck1.metadata,...
                JudgementsInfoToCheck2.matrix,...
                JudgementsInfoToCheck2.metadata,Rule3Ei);
        catch RepetitionsE
            Rule3E = addCause(Rule3E,RepetitionsE);
        end
    end
end

%% THROW THE MASTER EXCEPTION IF THE RULE WAS BROKEN
% If the master exception contains a cause, throw it out of this function.
if ~isempty(Rule3E.cause)
    throw(Rule3E)
end

end

%% FUNCTION TO CHECK JUDGEMENTS AGAINST IMPOSE JUDGEMENT RULE 4
function check_judgement_rule_4(~,JudgementsInclude,JudgementsInfo)
% This helper validates judgements againts rule 4.
% Rule 4 is that raw observable fixes must be continuous (with no gaps) and
% specified from the start of the forecast horizon.
%
% INPUTS:
%   -> JudgementsInclude: structure of indicators for each judgement type
%   -> JudgementsInfo: structure of information about judgements being made
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> overlay_one_incomplete_dataset_on_to_another
%   -> generate_MAPS_exception_and_add_as_cause

%% SETUP A MASTER EXCEPTION
% Create a master exception for rule 1 to add exceptions to as encountered
% below.
masterErrId = ['MAPS:',mfilename,':Rule4ValidationFailure'];
Rule4E = generate_MAPS_exception(masterErrId);

%% COMPUTE A COMBINED SET OF RAW OBSERVABLE FIXES
% Compute a combined set of raw observable fixes. If there exists only one
% type of raw observable fix, then the combined set is that type.
if JudgementsInclude.anticipatedRawObservableFixes && ...
        JudgementsInclude.unanticipatedRawObservableFixes
    [YtildejMat,YtildejMetadata] = ...
        overlay_one_incomplete_dataset_on_to_another(...
        JudgementsInfo.anticipatedRawObservableFixes.matrix,...
        JudgementsInfo.anticipatedRawObservableFixes.metadata,...
        JudgementsInfo.unanticipatedRawObservableFixes.matrix,...
        JudgementsInfo.unanticipatedRawObservableFixes.metadata);
elseif JudgementsInclude.anticipatedRawObservableFixes
    YtildejMat = JudgementsInfo.anticipatedRawObservableFixes.matrix;
    YtildejMetadata = ...
        JudgementsInfo.anticipatedRawObservableFixes.metadata;
elseif JudgementsInclude.unanticipatedRawObservableFixes
    YtildejMat = JudgementsInfo.unanticipatedRawObservableFixes.matrix;
    YtildejMetadata = ...
         JudgementsInfo.unanticipatedRawObservableFixes.metadata;
end

%% CHECK FOR VIOLATIONS OF RULE 4
% If any raw observable judgements exist in the judgement input, then
% compute the index numbers of any fixes that either do not begin in the
% first period of fix or that are discontinuous. Add causes to the master
% exception if either of these are non-empty.
if JudgementsInclude.anticipatedRawObservableFixes || ...
        JudgementsInclude.unanticipatedRawObservableFixes
    YtildejNaNinds = isnan(YtildejMat);
    BadInds.YtildeNotFromStartInds = find(YtildejNaNinds(:,1));
    BadInds.YtildeDiscontinuousInds = find(...
        any(diff(YtildejNaNinds,1,2)==-1,2));
    badIndTypes = {'NotFromStart';'Discontinuous'};
    nBadIndTypes = size(badIndTypes,1);
    for iType = 1:nBadIndTypes
        iBadIndType = ['Ytilde',badIndTypes{iType},'Inds'];
        iBadInds = BadInds.(iBadIndType);
        if ~isempty(iBadInds)
            errId = [masterErrId,':',badIndTypes{iType}];
            BadIndE = generate_MAPS_exception(errId);          
            niBadInds = size(iBadInds,1);
            for iiBadInd = 1:niBadInds
                BadIndE = generate_MAPS_exception_and_add_as_cause(...
                    BadIndE,[masterErrId,':Instance'],...
                    YtildejMetadata(iBadInds(iiBadInd)));
            end
            Rule4E = addCause(Rule4E,BadIndE);
        end
    end
end

%% THROW THE MASTER EXCEPTION IF THE RULE WAS BROKEN
% If the master exception contains a cause, throw it out of this function.
if ~isempty(Rule4E.cause)
    throw (Rule4E)
end

end

%% FUNCTION TO CHECK JUDGEMENTS AGAINST IMPOSE JUDGEMENT RULE 7
function check_judgement_rule_5(~,JudgementsInclude,JudgementsInfo)
% This helper validates judgements againts rule 5.
% Rule 5 is that the same shock cannot have judgement applied to it and be 
% used to implement a fix in the same period.
%
% INPUTS:   
%   -> JudgementsInclude: structure of indicators for each judgement type
%   -> JudgementsInfo: structure of information about judgements being made
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> validate_judgement_datasets_do_not_overlap (sub-function)

%% SETUP A MASTER EXCEPTION
% Create a master exception for rule 5 to add exceptions to as encountered
% below.
masterErrId = ['MAPS:',mfilename,':Rule5ValidationFailure'];
Rule5E = generate_MAPS_exception(masterErrId);

%% DEFINE JUDGEMENT PAIRS TO CHECK
% Define all combinations of shock judgements to check for overlaps.
judgementPairsToCheck = {
    'anticipatedShockJudgements'    'shocksForAnticipatedFixes'
    'unanticipatedShockJudgements'  'shocksForUnanticipatedFixes'};

%% RUN-THROUGH JUDGEMENTS TO CHECK
% If both judgement types are included in the judgements input, then check
% for repetitions. If any are found then add a cause to the master
% exception using the sub-function below.
nJudgementPairsToCheck = size(judgementPairsToCheck,1);
for iPair = 1:nJudgementPairsToCheck
    judgementToCheck1 = judgementPairsToCheck{iPair,1};
    judgementToCheck2 = judgementPairsToCheck{iPair,2};
    if JudgementsInclude.(judgementToCheck1) ...
            && JudgementsInclude.(judgementToCheck2)
        try
            errId = [masterErrId,':Instance'];
            Rule5Ei = generate_MAPS_exception(...
                errId,{judgementToCheck1 judgementToCheck2});
            validate_judgement_datasets_do_not_overlap(...
                JudgementsInfo.(judgementToCheck1).matrix,...
                JudgementsInfo.(judgementToCheck1).metadata,...
                JudgementsInfo.(judgementToCheck2).matrix,...
                JudgementsInfo.(judgementToCheck2).metadata,Rule5Ei);
        catch RepetitionsE
            Rule5E = addCause(Rule5E,RepetitionsE);
        end
    end
end

%% THROW THE MASTER EXCEPTION IF THE RULE WAS BROKEN
% If the master exception contains a cause, throw it out of this function.
if ~isempty(Rule5E.cause)
    throw(Rule5E)   
end

end

%% HELPER TO VALIDATE THAT TWO JUDGEMENTS DATASETS ARE NON-OVERLAPPING
function validate_judgement_datasets_do_not_overlap(...
    j1Mat,j1Metadata,j2Mat,j2Metadata,RepetitionE)
% This helper finds any overlap in the judgement datasets input.
% It throws an exception detailing the variables that are overlapping and
% the periods (forecast quarters) in which they overlap.
%
% INPUTS:   
%   -> j1Mat: first judgement dataset
%   -> j1Metadata: model metadata for the first judgement dataset
%   -> j2Mat: second judgement dataset
%   -> j2Metadata: model metadata for the second judgement dataset
%   -> RepetitionE: exception to add repetition causes to
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause    

%% FIND OVERLAPPING MODEL METADATA
% Use the MATLAB intersect function to find the model metadata that
% overlaps in the two judgements datasets.
[repMetadata,j1RepMetadataInds,j2RepMetadataInds] = intersect(...
    j1Metadata,j2Metadata);

%% COMPUTE PERIODS OF ANY OVERLAPS & CREATE EXCEPTIONS
% For each of the overlapping variables compute the periods of any
% overlaps. If there are any overlapping periods add a cause to a master
% exception detailing the variables that overlap and the periods (forecast 
% quarters) where those overalaps occurred.
if ~isempty(repMetadata)
    nReps = size(repMetadata,1);
    for iRep = 1:nReps
        iRepPeriodsOfReps = find(...
            ~isnan(j1Mat(j1RepMetadataInds(iRep),:))...
            .*~isnan(j2Mat(j2RepMetadataInds(iRep),:)));
        if ~isempty(iRepPeriodsOfReps)
            errId = ['MAPS:',mfilename,':InvalidOverlap'];
            errArgs = cell(1,2);
            errArgs{1} = j1Metadata{j1RepMetadataInds(iRep)};
            errArgs{2} = num2str(iRepPeriodsOfReps);
            RepetitionE = generate_MAPS_exception_and_add_as_cause(...
                RepetitionE,errId,errArgs);
        end
    end
end

%% THROW THE MASTER EXCEPTION IF ANY REPETITIONS WERE FOUND
% If the master exception contains a cause, throw it out of this function.
if ~isempty(RepetitionE.cause)
    throw(RepetitionE);
end

end

%% FUNCTION TO CREATE JUDGEMENT INFO
function JudgementsInfo = create_judgements_info_structure(...
    Model,Judgements,JudgementsInclude)
% This helper creates a judgement structure with all judgement information.
% It translates the judgements input into a new structure that is operated
% on by other sub-functions in the macro.
%
% INPUTS:
%   -> Model: MAPS LSS model structure
%   -> Judgements: judgements dataset structure
%   -> JudgementsInclude: structure describing the judgements content
%       - anticipatedModelVariableFixes: true or false
%       - anticipatedModelObservableFixes: true or false
%       - anticipatedRawObservableFixes: true or false
%       - shocksForAnticipatedFixes: true or false
%       - unanticipatedModelVariableFixes: true or false
%       - unanticipatedModelObservableFixes: true or false
%       - unanticipatedRawObservableFixes: true or false
%       - shocksForUnanticipatedFixes: true or false
%       - anticipatedShockJudgements: true or false
%       - unanticipatedShockJudgements true or false
%       - timeVaryingTrendJudgements: true or false
%
% OUTPUTS:
%   -> JudgementsInfo: new structure with as many fields as indicated true
%      in the judgement content descriptor input, each of which contain:
%       - metadata: judgement metadata (from judgements dataset)
%       - cell: judgement data as a cell array (from judgements dataset)
%       - modelMetadata: complete model metadata for that variable type
%       - indices: model index numbers equivalent to the metadata
%       - matrix: matrix equivalent to the cell vectors
%
% CALLS:
%   -> get_judgement_info_config (sub-function)
%   -> unpack_model
%   -> convert_zeros_to_NaNs (sub-function)
%   -> generate_MAPS_exception
%   -> lookup_model_index_numbers
%   -> generate_MAPS_exception_and_add_cause
%   -> convert_cell_array_of_vectors_to_matrix_equivalent

%% COMPLETE BASIC COMPONENTS OF THE STRUCTURE
% Call a sub-function which defines some configuration for the structure,
% including the name to associate that type of judgement with (as it
% appears in the judgements descriptor structure), the association with the
% input judgement structure and the model metadata to associate with that 
% type.
judgementInfoConfig = get_judgement_info_config;
allJudgementFields = judgementInfoConfig(:,1);
nAllFields = size(allJudgementFields,1);
for iField = 1:nAllFields
    if JudgementsInclude.(allJudgementFields{iField})
        JudgementsInfo.(allJudgementFields{iField}).metadata = ...
            judgementInfoConfig{iField,2}(Judgements);
        JudgementsInfo.(allJudgementFields{iField}).cell = ...
            judgementInfoConfig{iField,3}(Judgements);
        JudgementsInfo.(allJudgementFields{iField}).modelMetadata = ...
            unpack_model(Model,judgementInfoConfig(iField,4));
    end
end

%% CONVERT ZEROS TO NANS IN THE SHOCK USAGE VECTORS
% Convert the zero indicators for shock usage to NaNs so that they are
% consistent with and can be operated on in a similar to the judgements
% vectors (where non-usage is represented by a NaN).
shocksForFixesFields = {'shocksForAnticipatedFixes';
    'shocksForUnanticipatedFixes'};
nShocksForFixesFields = size(shocksForFixesFields,1);
for iField = 1:nShocksForFixesFields
    if JudgementsInclude.(shocksForFixesFields{iField})
        JudgementsInfo.(shocksForFixesFields{iField}).cell = ...
            convert_zeros_to_NaNs(...
            JudgementsInfo.(shocksForFixesFields{iField}).cell);
    end
end

%% CONVERT METADATA TO INDEX NUMBERS
% Convert the model metadata to index numbers using the model index lookup 
% function. Add any exceptions encountered as causes to a master exception,
% which is thrown at the end if non-empty. As this loops through, compute
% the dimensions (in number of periods) of the judgements, storing the
% maximum dimension encountered at each iteration.
masterErrId = ['MAPS:',mfilename,':InvalidJudgementMetadata'];
JmetadataE = generate_MAPS_exception(masterErrId);
judgementFields = fieldnames(JudgementsInfo);
nFields = size(judgementFields,1);
S = 0;
for iField = 1:nFields
    try
        JudgementsInfo.(judgementFields{iField}).indices = ...
            lookup_model_index_numbers(...
            JudgementsInfo.(judgementFields{iField}).modelMetadata,...
            JudgementsInfo.(judgementFields{iField}).metadata);
        Sj = max(cellfun('size',...
            JudgementsInfo.(judgementFields{iField}).cell,2));
        S = max(Sj,S);
    catch ModelIndexE
        errId = [masterErrId,':Instance'];
        JmetadataEi = generate_MAPS_exception_and_add_cause(ModelIndexE,...
            errId,allJudgementFields(iField));
        JmetadataE = addCause(JmetadataE,JmetadataEi);
    end
end
if ~isempty(JmetadataE.cause)
    throw(JmetadataE)
end

%% CONVERT THE CELL ARRAYS TO MATRICES
% Convert the cell arrays of judgement values to matrices of homegenous
% colum dimension (as given by the maximum judgement length calculated
% above).
for iField = 1:nFields
    JudgementsInfo.(judgementFields{iField}).matrix = ...
        convert_cell_array_of_vectors_to_matrix_equivalent(...
        JudgementsInfo.(judgementFields{iField}).cell,S);
end    
    

end

%% HELPER FUNCTION TO UNPACK CELL INPUTS
function judgementInfoConfig = get_judgement_info_config
% This helper contains a config for construction of judgement info.
% It associates several pieces of information with a judgement type as a
% helper for constructing the judgement info structure which is used to aid
% validation of the judgement against the impose judgement rules.
%
% INPUTS:   
%   -> none
%
% OUTPUTS:  
%   -> judgementInfoConfig: cell array of configuration information
%
% CALLS:
%   -> none

%% DEFINE THE CELL ARRAY CONFIGURATION
% Define the cell array configuration with each separate judgement type on
% a different row. Column 1 represents the judgement type; column 2 the
% location of the judgement metadata in the input; column 3 the location of
% the judgement cell array of vectors in the input; column 4 associates a 
% judgement with a particular type of model metadata.
judgementInfoConfig = {...
    'anticipatedModelVariableFixes' ...
        @(Judgements) Judgements.AnticipatedFixes.modelVariables(:,1) ...
        @(Judgements) Judgements.AnticipatedFixes.modelVariables(:,2) ...
        'xMnems';   
    'anticipatedModelObservableFixes' ...
        @(Judgements) Judgements.AnticipatedFixes.modelObservables(:,1) ...
        @(Judgements) Judgements.AnticipatedFixes.modelObservables(:,2) ...
        'Ymnems';
    'anticipatedRawObservableFixes' ...
        @(Judgements) Judgements.AnticipatedFixes.rawObservables(:,1) ...
        @(Judgements) Judgements.AnticipatedFixes.rawObservables(:,2) ...
        'YtildeMnems';
    'shocksForAnticipatedFixes' ...
        @(Judgements) Judgements.AnticipatedFixes.shockUsages(:,1) ...
        @(Judgements) Judgements.AnticipatedFixes.shockUsages(:,2) ...
        'zMnems';
    'unanticipatedModelVariableFixes' ...
        @(Judgements) Judgements.UnanticipatedFixes.modelVariables(:,1) ...
        @(Judgements) Judgements.UnanticipatedFixes.modelVariables(:,2) ...
        'xMnems';
    'unanticipatedModelObservableFixes' ...
        @(Judgements) Judgements.UnanticipatedFixes.modelObservables(:,1) ...
        @(Judgements) Judgements.UnanticipatedFixes.modelObservables(:,2) ...
        'Ymnems';
    'unanticipatedRawObservableFixes' ...
        @(Judgements) Judgements.UnanticipatedFixes.rawObservables(:,1) ...
        @(Judgements) Judgements.UnanticipatedFixes.rawObservables(:,2) ...
        'YtildeMnems';
    'shocksForUnanticipatedFixes' ...
        @(Judgements) Judgements.UnanticipatedFixes.shockUsages(:,1) ...       
        @(Judgements) Judgements.UnanticipatedFixes.shockUsages(:,2) ...
        'zMnems';
    'anticipatedShockJudgements' ...
        @(Judgements) Judgements.Shocks.anticipated(:,1) ...
        @(Judgements) Judgements.Shocks.anticipated(:,2) ...
        'zMnems';
    'unanticipatedShockJudgements' ...
        @(Judgements) Judgements.Shocks.unanticipated(:,1) ...
        @(Judgements) Judgements.Shocks.unanticipated(:,2) ...
        'zMnems'
    'timeVaryingTrendJudgements' ...
        @(Judgements) Judgements.timeVaryingTrends(:,1) ...
        @(Judgements) Judgements.timeVaryingTrends(:,2) ...
        'etatMnems'};
    
end

%% FUNCTION TO CONVERT THE ZEROS IN FIX SHOCK INDICATORS TO NANS
function ziCellWithNaNs = convert_zeros_to_NaNs(ziCellWithZeros)
% This helper converts all NaNs in a cell array of vectors to 0.
% It is used to make the indicator data for shock usage in the 
% implementation of fixes mimic the other types of judgement data where NaN
% signifies no usage. This allows the same functions to operate on all
% types of judgement data.
%
% INPUTS:   
%   -> ziCell: cell array of vectors with shock usage 0,1 indicators
%
% OUTPUTS:  
%   -> ziCell: cell array of vectors with shock usage NaN,1 indicators
%
% CALLS:
%   -> none

%% REPLACE THE ZEROS WITH NANS
% Loop over each of the numeric row vectors, replacing the zeros with NaNs
% in each vector.
ziCellWithNaNs = ziCellWithZeros;
nzi = size(ziCellWithZeros,1);
for izi = 1:nzi
    ziCellWithNaNs{izi}(ziCellWithZeros{izi}==0) = NaN;
end

end

%% FUNCTION TO VALIDATE JUDGEMENTS DATA STRUCTURE INPUT
function JudgementsInclude = interrogate_judgements_structure(...
    Model,Judgements)
% This helper validates the judgements input & determines its content. 
% It throws an exception if the inputs are invalid or inconsistent with the
% model input. Otherwise, it constructs a structure describing which types 
% of judgement are included in the input.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> Judgements: judgements dataset structure
%
% OUTPUTS:  
%   -> JudgementsInclude: structure describing the judgements content
%       - anticipatedModelVariableFixes: true or false
%       - anticipatedModelObservableFixes: true or false
%       - anticipatedRawObservableFixes: true or false
%       - shocksForAnticipatedFixes: true or false
%       - unanticipatedModelVariableFixes: true or false
%       - unanticipatedModelObservableFixes: true or false
%       - unanticipatedRawObservableFixes: true or false
%       - shocksForUnanticipatedFixes: true or false
%       - anticipatedShockJudgements: true or false
%       - unanticipatedShockJudgements true or false
%       - timeVaryingTrendJudgements: true or false
%
% CALLS:
%   -> generate_MAPS_exception
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause
%   -> generate_MAPS_exception_and_add_cause
%   -> validate_fixes_structure_input (sub-function)
%   -> validate_judgement_cell_array (sub-function)
%   -> unpack_model
     
%% SETUP INDICATOR VARIABLES FOR THE EXISTENCE OF EACH JUDGEMENT TYPE
% Setup an indicator in the output structure for the existence of all the
% possible fields in the original judgements input.
JudgementsInclude.anticipatedModelVariableFixes = false;
JudgementsInclude.anticipatedModelObservableFixes = false;
JudgementsInclude.anticipatedRawObservableFixes = false;
JudgementsInclude.shocksForAnticipatedFixes = false;
JudgementsInclude.unanticipatedModelVariableFixes = false;
JudgementsInclude.unanticipatedModelObservableFixes = false;
JudgementsInclude.unanticipatedRawObservableFixes = false;
JudgementsInclude.shocksForUnanticipatedFixes = false;
JudgementsInclude.anticipatedShockJudgements = false;
JudgementsInclude.unanticipatedShockJudgements = false;
JudgementsInclude.timeVaryingTrendJudgements = false;

%% CHECK BASIC SHAPE OF JUDGEMENTS INPUT
% Throw an exception if the judgements input does not contain any of the
% expected fields.
if ~isfield(Judgements,'AnticipatedFixes') && ...
        ~isfield(Judgements,'UnanticipatedFixes') && ...
        ~isfield(Judgements,'Shocks') && ...
        ~isfield(Judgements,'timeVaryingTrends')
    errId = ['MAPS:',mfilename,':MissingField'];
    generate_and_throw_MAPS_exception(errId);
end

%% CREATE EXCEPTION TO ADD JUDGEMENT INPUT EXCEPTIONS ENCOUNTERED BELOW
% Setup a master exception to add any causes as encountered below.
masterErrId = ['MAPS:',mfilename,':BadJudgementsInput'];
JinputE = generate_MAPS_exception(masterErrId);

%% CHECK FIXES INPUT
% Check for the validity of the fixes input (if it exists). The fixes
% input is invalid if: it is not a structure; it does not contain either of
% the expected fix type fields (Anticipated or Uannaticipated); or if the
% fields in these fix types are invalid (as checked in the
% check_fixes_input sub-function). If a particular type of fix is valid,
% then the field in the indicator structure is set to true.
varTypes = {'ModelVariable';'ModelObservable';'RawObservable'};
nVarTypes = size(varTypes,1);
fixTypes = {'Anticipated';'Unanticipated'};
nfixTypes = size(fixTypes,1);
for iType = 1:nfixTypes
    iFixType = [fixTypes{iType},'Fixes'];
    if isfield(Judgements,iFixType)
        try
            [existModelVariableFixes,existModelObservableFixes,...
                existRawObservableFixes] = ...
                validate_fixes_structure_input(Judgements.(iFixType));      %#ok<NASGU>
            for iiType = 1:nVarTypes
                iiFixType = lower(fixTypes{iType});
                iiVarType = varTypes{iiType};
                iVarTypeFieldName = [iiFixType,iiVarType,'Fixes'];
                iVarTypeField = eval(['exist',iiVarType,'Fixes']);
                JudgementsInclude.(iVarTypeFieldName) = iVarTypeField;
            end
            iShocksFieldName = ['shocksFor',fixTypes{iType},'Fixes'];
            JudgementsInclude.(iShocksFieldName) = true;
        catch FixesE
            errId = [masterErrId,':BadField'];
            JinputEi = generate_MAPS_exception_and_add_cause(...
                FixesE,errId,{iFixType});
            JinputE = addCause(JinputE,JinputEi);
        end
    end
end

%% CHECK SHOCKS JUDGEMENT INPUT
% Check for the validity of the shock judgements input (if it exists). The
% input is invalid if: it is not a structure; it does not contain either of
% the expected fix type fields (anticipated or uannaticipated); or if the
% fields in these shock judgement types are invalid (as checked in the
% check_judgement_cell_array sub-function). If a particular shock 
% judgement type is valid, then the field in the indicator structure is set 
% to true.
if isfield(Judgements,'Shocks')
    if ~isstruct(Judgements.Shocks)
        errId = [masterErrId,':BadShocksField'];
        JinputE = generate_MAPS_exception_and_add_as_cause(JinputE,errId);
    elseif ~isfield(Judgements.Shocks,'anticipated') && ...
            ~isfield(Judgements.Shocks,'unanticipated')
        errId = [masterErrId,':MissingShocksField'];
        JinputE = generate_and_add_MAPS_exception_cause(JinputE,errId);
    else
        shockJudgementTypes = {'anticipated';'unanticipated'};
        nShockJudgementTypes = size(shockJudgementTypes,1);
        for iType = 1:nShockJudgementTypes
            iShockJudgementType = shockJudgementTypes{iType};
            if isfield(Judgements.Shocks,iShockJudgementType)
                try
                    validate_judgement_cell_array(...
                        Judgements.Shocks.(iShockJudgementType));
                    iShockJudgementFieldName = ...
                        [iShockJudgementType,'ShockJudgements'];
                    JudgementsInclude.(iShockJudgementFieldName) = true;
                catch JcellE
                    errId = [masterErrId,':BadField'];
                    JinputEi = generate_MAPS_exception_and_add_cause(...
                        JcellE,errId,{['Shocks.',iShockJudgementType]});
                    JinputE = addCause(JinputE,JinputEi);
                end
            end
        end
    end
end

%% CHECK TIME-VARYING TREND JUDGEMENTS INPUT
% Check for the validity of the time-varying trend judgements input (if it 
% exists). The input is invalid if: it is not a structure; or if the
% field is invalid.
if isfield(Judgements,'timeVaryingTrends')
    try
        validate_judgement_cell_array(Judgements.timeVaryingTrends);
        JudgementsInclude.timeVaryingTrendJudgements = true;
    catch JcellE
        errId = [masterErrId,':BadField'];
        JinputEi = generate_MAPS_exception_and_add_cause(...
            JcellE,errId,{'timeVaryingTrends'});
        JinputE = addCause(JinputE,JinputEi);
    end
end

%% VALIDATE JUDGEMENTS AGAINST MODEL INPUTS
% Unpack the relevant model characteristics fields. Throw an error if: the
% judgement input includes anticipated judgements but the model is backward
% looking; the judgement input includes judgements in model observable
% space but the input model does not have measurement equations; the
% judgement input includes judgements in raw observable space but the model
% does not have data transformation equations; the judgement input includes
% time-varying deterministic trend judgements but the input model does not 
% have time-varying trends in its data transformation equations.
[modelIsForwardLooking,modelHasMeasurementEqs,...
    modelHasDataTransformationEqs,modelHasTimeVaryingTrends] = ...
    unpack_model(...
    Model,{'modelIsForwardLooking','modelHasMeasurementEqs',...
    'modelHasDataTransformationEqs','modelHasTimeVaryingTrends'});
if JudgementsInclude.anticipatedModelVariableFixes || ...
        JudgementsInclude.anticipatedModelObservableFixes || ...
        JudgementsInclude.anticipatedRawObservableFixes || ...
        JudgementsInclude.anticipatedShockJudgements
    if ~modelIsForwardLooking
        errId = [masterErrId,':UnexpectedAnticipatedJudgement'];
        JinputE = generate_MAPS_exception_and_add_as_cause(JinputE,errId);
    end
end
if JudgementsInclude.anticipatedModelObservableFixes || ...
        JudgementsInclude.unanticipatedModelObservableFixes
   if ~modelHasMeasurementEqs
       errId = [masterErrId,':UnexpectedModelObservableJudgement'];
       JinputE = generate_MAPS_exception_and_add_as_cause(JinputE,errId);
   end
end
if JudgementsInclude.anticipatedRawObservableFixes || ...
        JudgementsInclude.unanticipatedRawObservableFixes
    if ~modelHasDataTransformationEqs
        errId = [masterErrId,':UnexpectedRawObservableJudgement'];
        JinputE = generate_MAPS_exception_and_add_as_cause(JinputE,errId);
    end
end
if JudgementsInclude.timeVaryingTrendJudgements 
    if ~modelHasTimeVaryingTrends
        errId = [masterErrId,':UnexpectedTimeVaryingTrendJudgement'];
        JinputE = generate_MAPS_exception_and_add_as_cause(JinputE,errId);
    end
end

%% THROW ANY EXCEPTIONS ENCOUNTERED
% If the master exception contains any causes, then it is thrown and the
% function is terminated.
if ~isempty(JinputE.cause)
    throw(JinputE)
end

end

%% HELPER FUNCTION TO CHECK FIXES INPUT
function [existModelVariableFixes,existModelObservableFixes,...
    existRawObservableFixes] = validate_fixes_structure_input(Fixes)
% check_fixes_input validates the content of fields in the fixes input. 
% It throws ane exception if the input is not as excpeted and returns
% information about the content of the fixes otherwise.
%
% INPUTS:   
%   -> Fixes: one of the fixes components of the judgement input
%
% OUTPUTS:  
%   -> existModelVariableFixes: true if there are any fixes specified in
%      model variable space, false otherwise
%   -> existModelObservableFixes: true if there are any fixes specified in
%      model observable space, false otherwise
%   -> existRawObservableFixes: true if there are any fixes specified in
%      raw observable space, false otherwise
%
% CALLS:
%   -> generate_MAPS_exception
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception_and_add_cause
%   -> generate_MAPS_exception_and_add_as_cause
%   -> validate_judgement_cell_array (sub-function)
%   -> validate_shock_usage_cell_array (sub-function)

%% CHECK THE SHAPE OF THE INPUT
% Check that the component of the fixes input is a structure. If not throw
% an exception.
if ~isstruct(Fixes)
    errId = ['MAPS:',mfilename,':BadFixesInput'];
    generate_and_throw_MAPS_exception(errId);
end

%% CREATE MASTER EXCEPTION
% Create a master exception to add any causes encountered below.
masterErrId = ['MAPS:',mfilename,':BadFixesStructure'];
FixesE = generate_MAPS_exception(masterErrId);

%% CHECK SHOCK USAGE FOR THE FIXES
% Check that the fixes input have shocks associated with them. If not add
% an exception. If they do check the content of the shock usage in
% sub-functions defined below.
if ~isfield(Fixes,'shockUsages')
    errId = [masterErrId,':MissingShockUsagesField'];
    FixesE = generate_MAPS_exception_and_add_as_cause(FixesE,errId);
else
   try
       validate_judgement_cell_array(Fixes.shockUsages)
       validate_shock_usage_cell_array(Fixes.shockUsages);
   catch JcellE  
       errId = [masterErrId,':BadFieldInFixes'];
       FixesEi = generate_MAPS_exception_and_add_cause(...
           JcellE,errId,{'shockUsages'});
       FixesE = addCause(FixesE,FixesEi);
   end
end

%% CHECK THE FIXES
% Check that there are fixes specified in any of the spaces. If there are
% not add an exception. If yes, set the indicator outputs as appropriate
% and check that the content of the fixes in each space is as expected.
if ~isfield(Fixes,'modelVariables') && ...
        ~isfield(Fixes,'modelObservables') && ...
        ~isfield(Fixes,'rawObservables')
    errId = [masterErrId,':MissingFieldInFixes'];
    FixesE = generate_MAPS_exception_and_add_as_cause(FixesE,errId);
else
    varTypes = {'modelVariables';'modelObservables';'rawObservables'};
    nVarTypes = size(varTypes,1);
    existInInput = false(nVarTypes,1);
    for iType = 1:nVarTypes
        iVarType = varTypes{iType};
        if isfield(Fixes,iVarType)
            existInInput(iType) = true;
            try
                validate_judgement_cell_array(Fixes.(iVarType));
            catch JcellE
                errId = [masterErrId,':BadFieldInFixes'];
                FixesEi = generate_MAPS_exception_and_add_cause(...
                    JcellE,errId,{iVarType});
                FixesE = addCause(FixesE,FixesEi);
            end
        end
    end
    existModelVariableFixes = existInInput(1);
    existModelObservableFixes = existInInput(2);
    existRawObservableFixes = existInInput(3);
end

%% THROW THE EXCEPTION
% Throw the master exception if any exceptions were encountered above.
if ~isempty(FixesE.cause)
    throw(FixesE)
end

end
     
%% HELPER FUNCTION TO VALIDATE THE CONTENT OF A JUDGEMENT INPUT CELL
function validate_judgement_cell_array(jCell)
% This helper validates the content of a judgement cell input.
% Specifically, it checks that the data input is a two column cell array
% eith metadata (strings) in the first column and row vectors of data in
% the second column.
%
% INPUTS:   
%   -> jCell: data from a (terminal) field in the Judgements input
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause

%% CHECK SHAPE OF JUDGEMENT CELL INPUT
% Throw an exception if the judgement input is not a two-column cell array.
if ~iscell(jCell) || ndims(jCell)~=2 || size(jCell,2)~=2
    errId = ['MAPS:',mfilename,':BadJudgementField'];
    generate_and_throw_MAPS_exception(errId);
end

%% CREATE A MASTER EXCEPTION
% Create a master exception to add causes found below if any errors are 
% encountered.
masterErrId = ['MAPS:',mfilename,':BadJudgementCell'];
JcellE = generate_MAPS_exception(masterErrId);

%% CHECK JUDGEMENT METADATA
% Check that the judgement metadata is in the form of a non-empty cell
% string array. Add an exception as cause if not.
if ~iscellstr(jCell(:,1)) || any(cellfun(@isempty,jCell(:,1)))
    errId = [masterErrId,':BadMetadata'];
    JcellE = generate_MAPS_exception_and_add_as_cause(JcellE,errId);
end

%% CHECK JUDGEMENT VALUES
% Check that the judgement values constitute numeric, at least partly 
% finite, non-empty, two-dimensional row vectors. If they do not add an 
% exception as cause with the row numbers as arguments in that exception.
allInd = (1:size(jCell,1))';
badInd = find(cell2mat(cellfun(@(x) (~isnumeric(x)|ndims(x)~=2|...
    size(x,1)~=1|isempty(x)),jCell(:,2),'UniformOutput',false)));
goodInd = find(~ismember(allInd,badInd));
badGoodIndLogicals = cell2mat(cellfun(@(x) ~any(isfinite(x(:))),...
    jCell(goodInd,2),'UniformOutput',false));
badInd = sort([badInd;goodInd(badGoodIndLogicals)]);
if ~isempty(badInd)
    errId = [masterErrId,':BadData'];
    JcellE = generate_MAPS_exception_and_add_as_cause(...
        JcellE,errId,{num2str(badInd')});
end

%% THROW EXCEPTION
% Throw the master exception if either of the causes have been added.
if ~isempty(JcellE.cause)
    throw(JcellE)
end

end
    
%% HELPER FUNCTION TO VALIDATE SHOCKS TO BE USED IN FIXES ARE VALID
function validate_shock_usage_cell_array(fixShocks)
% This helper validates the input fix shock usage data. 
% Specifically, it checks that all of the data associated with shock usage
% in a particular type of fix is either 0 or 1.
%
% INPUTS:   
%   -> fixShocks: 2 column cell array with meteadata in the first column &
%      the shock usage data for each shock being used in the second column
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_and_throw_MAPS_exception

%% COMPUTE THE INDICES OF ANY SHOCK INDICATORS THAT BREAK THE RULES
% Use cellfun to find any elements in each vector of the cell that differ
% from 0 or 1. Record their index numbers.
badInd = find(cell2mat(cellfun(@(x) (any(x~=0&x~=1)|~any(x)),...
    fixShocks(:,2),'UniformOutput',false)));

%% WRITE & THROW ANY EXCEPTIONS ENOCUNTERED
% If any index numbers were found above throw and exception with the index 
% numbers as arguments.
if ~isempty(badInd)
    errId = ['MAPS:',mfilename,':BadShockUsageCell'];
    generate_and_throw_MAPS_exception(errId,{num2str(badInd')});
end 
        
end