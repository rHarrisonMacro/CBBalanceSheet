function ModelFileSyntaxE = check_MAPS_model_file_metadata_syntax(...
    ModelFileSyntaxE,FileContents,FileKeywords,FileLineNumbers,...
    metadataSyntaxChecksConfig)
% This function checks that metadata in a MAPS model file meets the rules. 
% It adds a cause to an exception for each separate rule that is broken.
% The rules relate to the validity of metdata in a MAPS model file from 
% uniqueness of variable and parameter mnemonics to the content of metadata 
% about the model.  
%
% INPUTS:
%   -> ModelFileSyntaxE: an exception to add cuases to
%   -> FileContents: structure containing all the information in the model
%      file deconstructed into its constituent parts
%   -> FileKeywords: structure (with the same fields as the above) 
%      containing the keywords used in the model file
%   -> FileLineNumbers: structure (with the same fields as the above) 
%      containing the line numbers in the file from which the information
%      was taken 
%   -> metadataSyntaxChecksConfig: cell array of configuration information
%      for the checks
%
% OUTPUTS:
%   -> ModelFileSyntaxE: updated exception
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> check_mnemonics_are_unique (sub-function)
%   -> check_mnemonics_are_used (sub-function)
%   -> check_mnemonics_are_continuous_expressions (sub-function)
%   -> check_mnemonic_string_lengths (sub-function)
%   -> check_mnemonics_are_not_banned (sub-function)
%   -> check_mnemonics_do_not_contain_banned_content (sub-function)
%   -> check_names_are_unique_within_fields (sub-function)
%   -> check_names_are_non_overlapping_across_fields (sub-function)
%   -> check_name_string_lengths (sub-function)
%   -> check_names_do_not_contain_banned_content (sub-function)
%   -> check_model_metadata_fields (sub-function)
%   -> check_parameter_values_are_valid (sub-function)
%
% DETAILS:
%   -> This helper contains a range of model file (and add-on file) syntax
%      checks that are common regardless of the type of model being parsed.
%   -> It is called by syntax checker functions that are specific to a
%      model type to avoid having to repeat checks across model types.
%   -> For example, it checks that the mnemonics across an entire model
%      file are unique.
%   -> Each individual check (and associated rule) is carried out in
%      separate sub-functions below.
%   -> The precise content of a check (and whether it happens at all) is
%      controlled by a configuration which is passed as input to this 
%      function and which will vary across model types.
%
% NOTES:
%   -> See <> for details of model file and add-on file syntax checking in
%      MAPS.
%
% This version: 06/06/2011
% Author(s): Alex Haberis & Matt Waldron

%% CHECK INPUTS
% Complete some basic checks on the input, including that the number of
% input arguments is correct.
if nargin < 5
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~strcmp(class(ModelFileSyntaxE),'MException')
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(FileContents)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(FileKeywords)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~isstruct(FileLineNumbers)
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);
elseif ~iscell(metadataSyntaxChecksConfig)
    errId = ['MAPS:',mfilename,':BadInput5'];
    generate_and_throw_MAPS_exception(errId);    
end

%% RUN THROUGH CHECKS
% Run through each of the checks dictated in the configuration where the
% configuration is an nChecks*1 cell array with the name of the check (as
% matched by a function name below) in the first column.
nChecksToRun = size(metadataSyntaxChecksConfig,1);
for iCheck = 1:nChecksToRun
    try
        iCheckFun = str2func(metadataSyntaxChecksConfig{iCheck,1});
        iCheckFun(FileContents,FileKeywords,FileLineNumbers,...
            metadataSyntaxChecksConfig{iCheck,2});
    catch SyntaxFunE
        ModelFileSyntaxE = addCause(ModelFileSyntaxE,SyntaxFunE);
    end
end

end

%% FUNCTION TO CHECK THAT ALL MODEL MNEMONICS ARE UNIQUE
function check_mnemonics_are_unique(...
    FileContents,~,FileLineNumbers,mnemFieldsToCheck)                       %#ok<DEFNU>
% This metadata syntax helper checks that all mnemonics are unique.
% This avoids errors and conflicts when two or more mnemonics in the model 
% are the same and they become confused.
%
% INPUTS:   
%   -> FileContents: structure containing all model info
%   -> FileLineNumbers: structure containing the line numbers in the file
%   -> mnemFieldsToCheck: names of the mnemonic fields to check
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> unpack_and_combine_file_contents (sub-function)
%   -> find_repeated_strings
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause

%% COLLECT ALL MNEMONICS IN FILE & LINE NUMBERS IN ONE VECTOR
% Use a hlper function to combine all the mnemonics specified in the 
% configuration input that exist in the file contents into one column cell
% string array. This function also returns the line numbers in the file on 
% which each of those mnemonics can be found.
[allMnems,allMnemsLineNumbers] = unpack_and_combine_file_contents(...
    mnemFieldsToCheck,FileContents,FileLineNumbers);

%% CHECK FOR REPEATED MNEMONICS
% Call a helper function to check for and return any repeated (non-unique) 
% mnemonics. Note that this is a case insensitive search.
repeatedMnems = find_repeated_strings(allMnems,true);
nRepeatedMnems = size(repeatedMnems,1);

%% CREATE & THROW AN EXCEPTION WITH REPEATED MNEMONICS
% If any repetitions were found, create an exception explaining the rule
% and then add each instance of repetition with the line numbers where
% those repetitions appear as exception causes to the master exception.
% Throw the master exception after all causes have been added.
if nRepeatedMnems > 0
    masterErrId = ['MAPS:',mfilename,':MnemsNotUnique'];
    MnemsNotUniqueE = generate_MAPS_exception(masterErrId);
    errId = [masterErrId,':Instance'];
    errArgs = cell(1,2);
    for iRep = 1:nRepeatedMnems
        iRepLogicals = strcmpi(repeatedMnems{iRep},allMnems);
        errArgs{1} = repeatedMnems{iRep};
        errArgs{2} = num2str(allMnemsLineNumbers(iRepLogicals)');
        MnemsNotUniqueE = generate_MAPS_exception_and_add_as_cause(...
            MnemsNotUniqueE,errId,errArgs);
    end
    throw(MnemsNotUniqueE);
end

end

%% FUNCTION TO CHECK MNEMONICS ARE USED
function check_mnemonics_are_used(...
    FileContents,FileKeywords,FileLineNumbers,mnemUsageCheckConfig)         %#ok<DEFNU>
% This metadata syntax helper checks that all mnemonics are used.
% This helps users to keep their model files tidy and up-to-date.
%
% INPUTS:   
%   -> FileContents: structure containing all model info
%   -> FileKeywords: structure containing the keywords used in the file 
%   -> FileLineNumbers: structure containing the line numbers in the file
%   -> mnemUsageCheckConfig: configuration for the check with fields
%      names of the mnemonics to check & equations to check them in
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> unpack_and_combine_file_contents (sub-function)
%   -> split_equation
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause

%% UNPACK CONFIGURATION INFO
% Unpack the list of all the mnemonics to check and the equations to check
% they're used in from the first and second elements of the input 
% configuration cell.
mnemFieldsToCheck = mnemUsageCheckConfig{1};
eqFieldsToCheckThemIn = mnemUsageCheckConfig{2};

%% COLLECT ALL MNEMONICS IN FILE, LINE NUMBERS & KEYWORDS IN ONE VECTOR
% Use a hlper function to combine all the mnemonics specified in the 
% configuration input that exist in the file contents into one column cell
% string array. This function also returns the line numbers in the file on 
% which each of those mnemonics can be found and the keywords under which
% they appear.
[allMnems,allMnemsLineNumbers,allMnemsKeywords] = ...
    unpack_and_combine_file_contents(...
    mnemFieldsToCheck,FileContents,FileLineNumbers,FileKeywords);

%% BAIL OUT IF THERE ARE NO MNEMONICS TO CHECK
% This scenario is never likely to occur for NLBL or LSS models, but is
% entirely feasible for a decomp add-on file
if isempty(allMnems) || isequal (allMnems,{''})
    return;
end

%% COLLECT ALL EQUATIONS IN ONE VECTOR
% Use the same helper to combine all the equations specified in the 
% configuration input that exist in the file contents into one column cell
% string array.
allEqs = unpack_and_combine_file_contents(...
    eqFieldsToCheckThemIn,FileContents);

%% SPLIT ALL EQUATION TERMS OUT & COLLECT IN ONE VECTOR
% Split all the equation terms out for each of the equation separately and
% then comine into one column cell string array.
nAllEqs = size(allEqs,1);
allEqsSplitTerms = cell(nAllEqs,1);
for iEq = 1:nAllEqs
    allEqsSplitTerms{iEq} = split_equation(allEqs{iEq});
end
allEqsTerms = [allEqsSplitTerms{:}]';

%% FIND UNUNSED MNEMONICS
% Find the index numbers of all the mnemonics that are not found among the
% complete cell string array of equation terms using the MATLAB ismember
% function.
ununsedMnemInds = find(~ismember(allMnems,allEqsTerms));
nUnusedMnems = size(ununsedMnemInds,1);

%% CREATE & THROW AN EXCEPTION WITH UNUSED MNEMONICS
% If any unused mnemonics were found, create an exception explaining the 
% rule and then add each instance of unused mnemonic with the line number 
% where the mnemonic appears and the model file keyword under which it
% appears as exception causes to the master exception. Throw the master 
% exception after all causes have been added.
if nUnusedMnems > 0
    masterErrId = ['MAPS:',mfilename,':MnemsNotUsed'];
    MnemsNotUsedE = generate_MAPS_exception(masterErrId);    
    errId = [masterErrId,':Instance'];
    errArgs = cell(1,3);
    for iMnem = 1:nUnusedMnems
        iUnusedMnemInd = ununsedMnemInds(iMnem);
        errArgs{1} = allMnems{iUnusedMnemInd};
        errArgs{2} = num2str(allMnemsLineNumbers(iUnusedMnemInd));
        errArgs{3} = allMnemsKeywords{iUnusedMnemInd};
        MnemsNotUsedE = generate_MAPS_exception_and_add_as_cause(...
            MnemsNotUsedE,errId,errArgs);
    end
    throw(MnemsNotUsedE);
end

end

%% FUNCTION TO CHECK MNEMONICS ARE CONTINUOUS EXPRESSIONS
function check_mnemonics_are_continuous_expressions(...
    FileContents,FileKeywords,FileLineNumbers,mnemFieldsToCheck)                 %#ok<DEFNU>
% This metadata syntax helper checks that all mnemonics are continuous.
% This avoids errors in symbolic model creation when a mnemonic is defined
% as two or more expressions.
%
% INPUTS:   
%   -> FileContents: structure containing all model info
%   -> FileKeywords: structure containing the keywords used in the file 
%   -> FileLineNumbers: structure containing the line numbers in the file
%   -> mnemFieldsToCheck: names of the mnemonic fields to check
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> find_fields_in_model (sub-function)

%% COLLECT ALL MNEMONICS IN FILE & LINE NUMBERS IN ONE VECTOR
% Use a hlper function to combine all the mnemonics specified in the 
% configuration input that exist in the file contents into one column cell
% string array. This function also returns the line numbers in the file on 
% which each of those mnemonics can be found.
[allMnems,allMnemsLineNumbers,allMnemsKeywords] = ...
    unpack_and_combine_file_contents(...
    mnemFieldsToCheck,FileContents,FileLineNumbers,FileKeywords);

%% CHECK FOR MNEMONICS WITH SPACE CHARACTERS
% Search all the mnemoncs for white space characters ('\s') using the
% MATLAB regexp command.
nonContinuousMnemInds = find(~cellfun(...
    @isempty,regexp(allMnems,'\s','match')));
nNonContinuousMnems = size(nonContinuousMnemInds,1);

%% CREATE & THROW AN EXCEPTION WITH UNUSED MNEMONICS
% If any unused mnemonics were found, create an exception explaining the 
% rule and then add each instance of unused mnemonic with the line number 
% where the mnemonic appears and the model file keyword under which it
% appears as exception causes to the master exception. Throw the master 
% exception after all causes have been added.
if nNonContinuousMnems > 0
    masterErrId = ['MAPS:',mfilename,':MnemsNotContinuous'];
    MnemsNotContinuousE = generate_MAPS_exception(masterErrId);    
    errId = [masterErrId,':Instance'];
    errArgs = cell(1,3);
    for iMnem = 1:nNonContinuousMnems
        iNonContinuousMnemInd = nonContinuousMnemInds(iMnem);
        errArgs{1} = allMnems{iNonContinuousMnemInd};
        errArgs{2} = num2str(allMnemsLineNumbers(iNonContinuousMnemInd));
        errArgs{3} = allMnemsKeywords{iNonContinuousMnemInd};
        MnemsNotContinuousE = generate_MAPS_exception_and_add_as_cause(...
            MnemsNotContinuousE,errId,errArgs);
    end
    throw(MnemsNotContinuousE);
end

end

%% FUNCTION TO CHECK THAT MNEMONIC STRINGS ARE NOT TOO LONG
function check_mnemonic_string_lengths(...
    FileContents,FileKeywords,FileLineNumbers,mnemLengthCheckConfig)        %#ok<DEFNU>
% This metadata syntax helper checks that mnemonic stings are not too long.
% This avoids errors in EASE and in publication of data from EASE.
%
% INPUTS:   
%   -> FileContents: structure containing all model info
%   -> FileKeywords: structure containing the keywords used in the file 
%   -> FileLineNumbers: structure containing the line numbers in the file
%   -> mnemLengthCheckConfig: configuration for the check with fields
%      names of the mnemonics to check & a maximum mnemonic length
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> unpack_and_combine_file_contents (sub-function)
%   -> generate_MAPS_exception 
%   -> generate_MAPS_exception_and_add_as_cause 

%% UNPACK CONFIGURATION INFO
% Unpack the list of all the mnemonics to check and the equations to check
% they're used in from the first and second elements of the input 
% configuration cell.
mnemFieldsToCheck = mnemLengthCheckConfig{1};
maxMnemLength = mnemLengthCheckConfig{2};

%% COLLECT ALL MNEMONICS IN FILE & LINE NUMBERS IN ONE VECTOR
% Use a hlper function to combine all the mnemonics specified in the 
% configuration input that exist in the file contents into one column cell
% string array. This function also returns the line numbers in the file on 
% which each of those mnemonics can be found.
[allMnems,allMnemsLineNumbers,allMnemsKeywords] = ...
    unpack_and_combine_file_contents(...
    mnemFieldsToCheck,FileContents,FileLineNumbers,FileKeywords);

%% CHECK FOR MNEMONICS WITH LENGTHS EXCEEDING THE MAXIMUM ALLOWED
% Find all the mnemonics which have a length exceeding the maximum length
% allowed (as covered by the configuration input).
tooLongMnemInds = find(cellfun('size',allMnems,2)>maxMnemLength);
nTooLongMnems = size(tooLongMnemInds,1);

%% CREATE & THROW AN EXCEPTION WITH INVALID MNEMONICS
% If any mnemonic strings with too many characters were found, create an 
% exception explaining the rule and then add each instance with the line 
% number where the mnemonic appears and the model file keyword under which 
% it appears as exception causes to the master exception. Throw the master 
% exception after all causes have been added.
if nTooLongMnems > 0
    masterErrId = ['MAPS:',mfilename,':MnemsTooLong'];
    masterErrArgs = {num2str(maxMnemLength)};
    MnemsTooLongE = generate_MAPS_exception(masterErrId,masterErrArgs);    
    errId = [masterErrId,':Instance'];
    errArgs = cell(1,3);
    for iMnem = 1:nTooLongMnems
        iMnemTooLongInd = tooLongMnemInds(iMnem);
        errArgs{1} = allMnems{iMnemTooLongInd};
        errArgs{2} = num2str(allMnemsLineNumbers(iMnemTooLongInd));
        errArgs{3} = allMnemsKeywords{iMnemTooLongInd};
        MnemsTooLongE = generate_MAPS_exception_and_add_as_cause(...
            MnemsTooLongE,errId,errArgs);
    end
    throw(MnemsTooLongE);
end

end

%% FUNCTION TO CHECK THAT MNEMONICS ARE NOT BANNED
function check_mnemonics_are_not_banned(...
    FileContents,FileKeywords,FileLineNumbers,mnemFieldsToCheck)            %#ok<DEFNU>
% This metadata syntax helper checks that mnemonics are not banned.
% This avoids errors in MAPS (and EASE) when mnemonics have the same names
% as mathematical operators, for example.
%
% INPUTS:   
%   -> FileContents: structure containing all model info
%   -> FileKeywords: structure containing the keywords used in the file 
%   -> FileLineNumbers: structure containing the line numbers in the file
%   -> mnemFieldsToCheck: names of the mnemonic fields to check
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> unpack_and_combine_file_contents (sub-function)
%   -> get_banned_mnemonics_config
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause

%% COLLECT ALL MNEMONICS IN FILE & LINE NUMBERS IN ONE VECTOR
% Use a hlper function to combine all the mnemonics specified in the 
% configuration input that exist in the file contents into one column cell
% string array. This function also returns the line numbers in the file on 
% which each of those mnemonics can be found.
[allMnems,allMnemsLineNumbers,allMnemsKeywords] = ...
    unpack_and_combine_file_contents(...
    mnemFieldsToCheck,FileContents,FileLineNumbers,FileKeywords);

%% GET BANNED MNEOMICS CONFIGURATION
% Get the list of banned mnemonics from the configuration file.
bannedMnems = get_banned_mnemonics_config;

%% CHECK FOR MNEMONICS THAT ARE BANNED
% Check that each mnemonic in the model is not banned and find the index
% numbers of those that are
bannedMnemInds = find(cell2mat(cellfun(...
    @(x) any(strcmpi(x,bannedMnems)),allMnems,'UniformOutput',false)));
nBannedMnems = size(bannedMnemInds,1);

%% CREATE & THROW AN EXCEPTION WITH BANNED MNEMONICS
% If any banned mnemonics were found, create an exception explaining the 
% rule and then add each instance with the line number where the mnemonic 
% appears and the model file keyword under which it appears as exception 
% causes to the master exception. Throw the master exception after all 
% causes have been added.
if nBannedMnems > 0
    masterErrId = ['MAPS:',mfilename,':BannedMnems'];
    masterErrArgs = bannedMnems;
    BannedMnemsE = generate_MAPS_exception(masterErrId,masterErrArgs);    
    errId = [masterErrId,':Instance'];
    errArgs = cell(1,3);
    for iMnem = 1:nBannedMnems
        iBannedMnemInd = bannedMnemInds(iMnem);
        errArgs{1} = allMnems{iBannedMnemInd};
        errArgs{2} = num2str(allMnemsLineNumbers(iBannedMnemInd));
        errArgs{3} = allMnemsKeywords{iBannedMnemInd};
        BannedMnemsE = generate_MAPS_exception_and_add_as_cause(...
            BannedMnemsE,errId,errArgs);
    end
    throw(BannedMnemsE);
end

end

%% FUNCTION TO CHECK THAT MNEMONICS DO NOT CONTAIN BANNED CONTENT
function check_mnemonics_do_not_contain_banned_content(...
    FileContents,FileKeywords,FileLineNumbers,mnemFieldsToCheck)            %#ok<DEFNU>
% This syntax helper checks that mnemonics do not contain banned content.
% This avoids errors in MAPS (and EASE) when mnemonics contain symbols that
% have special meaning like mathematical operators, for example.
%
% INPUTS:   
%   -> FileContents: structure containing all model info
%   -> FileKeywords: structure containing the keywords used in the file 
%   -> FileLineNumbers: structure containing the line numbers in the file
%   -> mnemFieldsToCheck: names of the mnemonic fields to check
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> unpack_and_combine_file_contents (sub-function)
%   -> get_banned_mnemonic_content_config
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause

%% COLLECT ALL MNEMONICS IN FILE & LINE NUMBERS IN ONE VECTOR
% Use a hlper function to combine all the mnemonics specified in the 
% configuration input that exist in the file contents into one column cell
% string array. This function also returns the line numbers in the file on 
% which each of those mnemonics can be found.
[allMnems,allMnemsLineNumbers,allMnemsKeywords] = ...
    unpack_and_combine_file_contents(...
    mnemFieldsToCheck,FileContents,FileLineNumbers,FileKeywords);

%% GET BANNED MNEOMICS CONTENT CONFIGURATION
% Get the list of banned mnemonic content from the configuration file.
bannedMnemContent = get_banned_mnemonic_content_config;

%% CHECK FOR MNEMONICS THAT CONTAIN BANNED CONTENT
% Compute logicals for each of the mnemonics describing whether or not each 
% banned symbol appears in the mnemonic in an 1*nBannedMnemContent cell 
% array of nAllMnems*1 logicals. Combine those together into one matrix.
% And then conpute the index numbers of the invalid mnemonics.
bannedMnemContentLogicalsCell = cellfun(@(x) ~cellfun(...
    @isempty,strfind(allMnems,x)),bannedMnemContent,...
    'UniformOutput',false);
bannedMnemContentLogicals = [bannedMnemContentLogicalsCell{:}];
mnemsWithBannedContentInds = find(any(bannedMnemContentLogicals,2));
nMnemsWithBannedContent = size(mnemsWithBannedContentInds,1);

%% CREATE & THROW AN EXCEPTION WITH MNEMONICS THAT CONTAIN BANNED CONTENT
% If any mnemonics with banned content were found, create an exception 
% explaining the rule and then add each instance with the line number where 
% the mnemonic appears and the model file keyword under which it appears as 
% exception causes to the master exception. Throw the master exception 
% after all causes have been added.
if nMnemsWithBannedContent > 0
    masterErrId = ['MAPS:',mfilename,':BannedMnemContent'];
    masterErrArgs = bannedMnemContent;
    BannedMnemContentE = generate_MAPS_exception(...
        masterErrId,masterErrArgs);    
    errId = [masterErrId,':Instance'];
    errArgs = cell(1,4);
    for iMnem = 1:nMnemsWithBannedContent
        iMnemWithBannedContentInd = mnemsWithBannedContentInds(iMnem);
        errArgs{1} = allMnems{iMnemWithBannedContentInd};
        errArgs{2} = num2str(...
            allMnemsLineNumbers(iMnemWithBannedContentInd));
        errArgs{3} = allMnemsKeywords{iMnemWithBannedContentInd};
        errArgs{4} = bannedMnemContent(...
            bannedMnemContentLogicals(iMnemWithBannedContentInd,:));
        BannedMnemContentE = generate_MAPS_exception_and_add_as_cause(...
            BannedMnemContentE,errId,errArgs);
    end
    throw(BannedMnemContentE);
end

end

%% FUNCTION TO CHECK THAT MODEL NAMES WITHIN A FIELD ARE UNIQUE
function check_names_are_unique_within_fields(...
    FileContents,FileKeywords,FileLineNumbers,nameFieldsToCheck)            %#ok<DEFNU>
% This metadata syntax helper checks that names within a field are unique.
% This avoids errors and conflicts when two or more names within a model 
% field are the same and they become confused (eg in equation 
% decompositions).
%
% INPUTS:   
%   -> FileContents: structure containing all model info 
%   -> FileKeywords: structure containing the keywords used in the file 
%   -> FileLineNumbers: structure containing the line numbers in the file
%   -> namesFieldsToCheck: names of the name fields to check
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> find_repeated_strings
%   -> generate_MAPS_exception_and_add_as_cause

%% SETUP A MASTER EXCEPTION
% Setup a master exception explaining the syntax rule to which causes can
% be added as appropriate below. Setup an identifier for the causes and an
% error arguments cell to add specific instances to.
masterErrId = ['MAPS:',mfilename,':NamesNotUniqueWithinField'];
NamesNotUniqueE = generate_MAPS_exception(masterErrId);
errId = [masterErrId,':Instance'];
errArgs = cell(1,3);

%% CHECK NAMES IN EACH FIELD
% For each field input, check that the field exists in the file contents.
% If it does, check that there are no repeated names (regardless of case).
% If there are, add causes to the master exception as appropriate.
nNameFields = size(nameFieldsToCheck,1);
for iField = 1:nNameFields
    iFieldName = nameFieldsToCheck{iField};
    if isfield(FileContents,iFieldName)
        namesToCheck = FileContents.(iFieldName);
        namesToCheckKeyword = FileKeywords.(iFieldName);
        namesToCheckLineNumbers = FileLineNumbers.(iFieldName);
        repeatedNames = find_repeated_strings(namesToCheck,true);
        nRepeatedNames = size(repeatedNames,1);
        if nRepeatedNames > 0                      
            for iRep = 1:nRepeatedNames
                iRepLogicals = strcmpi(repeatedNames{iRep},namesToCheck);
                errArgs{1} = repeatedNames{iRep};
                errArgs{2} = ...
                    num2str(namesToCheckLineNumbers(iRepLogicals)');
                errArgs{3} = namesToCheckKeyword;
                NamesNotUniqueE = ...
                    generate_MAPS_exception_and_add_as_cause(...
                    NamesNotUniqueE,errId,errArgs);
            end
        end
    end
end

%% THROW MASTER EXCEPTION
% If the master exception contains any causes, throw it out of the
% function.
if ~isempty(NamesNotUniqueE.cause)
    throw(NamesNotUniqueE);
end

end

%% FUNCTION TO CHECK THAT NAMES ARE NON OVERLAPPING ACROSS FIELDS
function check_names_are_non_overlapping_across_fields(...
    FileContents,FileKeywords,FileLineNumbers,namePairsToCheck)             %#ok<DEFNU>
% This metadata syntax helper checks that names across fields are unique.
% This avoids errors and conflicts when two or more names across model 
% field are the same and they become confused (eg in equation 
% decompositions).
%
% INPUTS:   
%   -> FileContents: structure containing all model info 
%   -> FileKeywords: structure containing the keywords used in the file 
%   -> FileLineNumbers: structure containing the line numbers in the file
%   -> namePairsToCheck: names of the fields to check
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause

%% SETUP A MASTER EXCEPTION
% Setup a master exception explaining the syntax rule to which causes can
% be added as appropriate below. Setup an identifier for the causes and an
% error arguments cell to add specific instances to.
masterErrId = ['MAPS:',mfilename,':NamesOverlapAcrossFields'];
NamesOverlapE = generate_MAPS_exception(masterErrId);
errId = [masterErrId,':Instance'];
errArgs = cell(1,5);

%% CHECK NAMES ACROSS FIELD PAIRS
% For each field pair input, check that both field exists in the file 
% contents. If they do, check that there are no repeated names. If there 
% are, add causes to the master exception as appropriate.
nNamePairs = size(namePairsToCheck,1);
for iPair = 1:nNamePairs
    iFieldName1 = namePairsToCheck{iPair,1};
    iFieldName2 = namePairsToCheck{iPair,2};
    if isfield(FileContents,iFieldName1) && ...
            isfield(FileContents,iFieldName2)
        namesToCheck1 = FileContents.(iFieldName1);
        namesToCheck2 = FileContents.(iFieldName2);
        namesToCheck1Keyword = FileKeywords.(iFieldName1);
        namesToCheck2Keyword = FileKeywords.(iFieldName2);
        namesToCheck1LineNumbers = FileLineNumbers.(iFieldName1);
        namesToCheck2LineNumbers = FileLineNumbers.(iFieldName2);
        iPairOverlapLogicals = ismember(...
            lower(namesToCheck1),lower(namesToCheck2));
        if any(iPairOverlapLogicals)
            iPairOverlapInds = find(iPairOverlapLogicals);
            niPairOverlapInds = size(iPairOverlapInds,1);
            for iInd = 1:niPairOverlapInds
                iOverlappingName = namesToCheck1{iPairOverlapInds(iInd)};               
                errArgs{1} = iOverlappingName;
                errArgs{2} = num2str(...
                    namesToCheck1LineNumbers(iPairOverlapInds(iInd)));
                errArgs{3} = namesToCheck1Keyword;
                errArgs{4} = num2str(namesToCheck2LineNumbers(...
                    strcmpi(iOverlappingName,namesToCheck2))');
                errArgs{5} = namesToCheck2Keyword;
                NamesOverlapE = ...
                    generate_MAPS_exception_and_add_as_cause(...
                    NamesOverlapE,errId,errArgs);
            end
        end
    end
end

%% THROW MASTER EXCEPTION
% If the master exception contains any causes, throw it out of the
% function.
if ~isempty(NamesOverlapE.cause)
    throw(NamesOverlapE);
end

end

%% FUNCTION TO CHECK THAT NAME STRINGS ARE NOT TOO LONG
function check_name_string_lengths(...
    FileContents,FileKeywords,FileLineNumbers,nameLengthCheckConfig)        %#ok<DEFNU>
% This metadata syntax helper checks that name stings are not too long.
% This avoids errors in EASE with storage of names in the database.
%
% INPUTS:   
%   -> FileContents: structure containing all model info
%   -> FileKeywords: structure containing the keywords used in the file 
%   -> FileLineNumbers: structure containing the line numbers in the file
%   -> nameLengthCheckConfig: configuration for the check with field
%      names of the names to check & a maximum name length
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> unpack_and_combine_file_contents (sub-function)
%   -> generate_MAPS_exception 
%   -> generate_MAPS_exception_and_add_as_cause 

%% UNPACK CONFIGURATION INFO
% Unpack the list of all the mnemonics to check and the equations to check
% they're used in from the first and second elements of the input 
% configuration cell.
nameFieldsToCheck = nameLengthCheckConfig{1};
maxNameLength = nameLengthCheckConfig{2};

%% COLLECT ALL NAMES IN FILE & LINE NUMBERS IN ONE VECTOR
% Use a helper function to combine all the names specified in the 
% configuration input that exist in the file contents into one column cell
% string array. This function also returns the line numbers in the file on 
% which each of those names can be found.
[allNames,allNamesLineNumbers,allNamesKeywords] = ...
    unpack_and_combine_file_contents(...
    nameFieldsToCheck,FileContents,FileLineNumbers,FileKeywords);

%% CHECK FOR NAMES WITH LENGTHS EXCEEDING THE MAXIMUM ALLOWED
% Find all the names which have a length exceeding the maximum length
% allowed (as covered by the configuration input).
tooLongNameInds = find(cellfun('size',allNames,2)>maxNameLength);
nTooLongNames = size(tooLongNameInds,1);

%% CREATE & THROW AN EXCEPTION WITH INVALID NAMES
% If any name strings with too many characters were found, create an 
% exception explaining the rule and then add each instance with the line 
% number where the name appears and the model file keyword under which 
% it appears as exception causes to the master exception. Throw the master 
% exception after all causes have been added.
if nTooLongNames > 0
    masterErrId = ['MAPS:',mfilename,':NamesTooLong'];
    masterErrArgs = {num2str(maxNameLength)};
    NamesTooLongE = generate_MAPS_exception(masterErrId,masterErrArgs);    
    errId = [masterErrId,':Instance'];
    errArgs = cell(1,3);
    for iName = 1:nTooLongNames
        iNameTooLongInd = tooLongNameInds(iName);
        errArgs{1} = allNames{iNameTooLongInd};
        errArgs{2} = num2str(allNamesLineNumbers(iNameTooLongInd));
        errArgs{3} = allNamesKeywords{iNameTooLongInd};
        NamesTooLongE = generate_MAPS_exception_and_add_as_cause(...
            NamesTooLongE,errId,errArgs);
    end
    throw(NamesTooLongE);
end

end

%% FUNCTION TO CHECK THAT NAMES DO NOT CONTAIN BANNED CONTENT
function check_names_do_not_contain_banned_content(...
    FileContents,FileKeywords,FileLineNumbers,nameFieldsToCheck)            %#ok<DEFNU>
% This syntax helper checks that names do not contain banned content.
% This avoids errors in MAPS (and EASE) when names contain symbols that
% have special meaning like mathematical operators, for example.
%
% INPUTS:   
%   -> FileContents: structure containing all model info
%   -> FileKeywords: structure containing the keywords used in the file 
%   -> FileLineNumbers: structure containing the line numbers in the file
%   -> nameFieldsToCheck: names of the name fields to check
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> unpack_and_combine_file_contents (sub-function)
%   -> get_banned_name_content_config
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause

%% COLLECT ALL NAMES IN FILE & LINE NUMBERS IN ONE VECTOR
% Use a hlper function to combine all the names specified in the 
% configuration input that exist in the file contents into one column cell
% string array. This function also returns the line numbers in the file on 
% which each of those names can be found and the keywords that apply.
[allNames,allNamesLineNumbers,allNamesKeywords] = ...
    unpack_and_combine_file_contents(...
    nameFieldsToCheck,FileContents,FileLineNumbers,FileKeywords);

%% GET BANNED NAME CONTENT CONFIGURATION
% Get the list of banned name content from the configuration file.
bannedNameContent = get_banned_name_content_config;

%% CHECK FOR NAMES THAT CONTAIN BANNED CONTENT
% For each banned symbol, compute logicals for each of the names describing
% whether or not each banned symbol appears in the mnemonic in an 
% 1*nBannedNameContent cell 
% array of nAllMnems*1 logicals. Combine those together into one matrix.
% And then conpute the index numbers of the invalid mnemonics.
bannedNameContentLogicalsCell = cellfun(@(x) ~cellfun(...
    @isempty,strfind(allNames,x)),bannedNameContent,...
    'UniformOutput',false);
bannedNameContentLogicals = [bannedNameContentLogicalsCell{:}];
namesWithBannedContentInds = find(any(bannedNameContentLogicals,2));
nNamesWithBannedContent = size(namesWithBannedContentInds,1);

%% CREATE & THROW AN EXCEPTION WITH NAMES THAT CONTAIN BANNED CONTENT
% If any names with banned content were found, create an exception 
% explaining the rule and then add each instance with the line number where 
% the name appears and the model file keyword under which it appears as 
% exception causes to the master exception. Throw the master exception 
% after all causes have been added.
if nNamesWithBannedContent > 0
    masterErrId = ['MAPS:',mfilename,':BannedNameContent'];
    masterErrArgs = bannedNameContent;
    BannedNameContentE = generate_MAPS_exception(...
        masterErrId,masterErrArgs);    
    errId = [masterErrId,':Instance'];
    errArgs = cell(1,4);
    for iName = 1:nNamesWithBannedContent
        iNameWithBannedContentInd = namesWithBannedContentInds(iName);
        errArgs{1} = allNames{iNameWithBannedContentInd};
        errArgs{2} = num2str(...
            allNamesLineNumbers(iNameWithBannedContentInd));
        errArgs{3} = allNamesKeywords{iNameWithBannedContentInd};
        errArgs{4} = bannedNameContent(...
            bannedNameContentLogicals(iNameWithBannedContentInd,:));
        BannedNameContentE = generate_MAPS_exception_and_add_as_cause(...
            BannedNameContentE,errId,errArgs);
    end
    throw(BannedNameContentE);
end

end

%% FUNCTION TO CHECK MODEL METADATA CONTENT
function check_model_metadata_fields(...
    FileContents,FileKeywords,FileLineNumbers,metadataFieldsCheckConfig)    %#ok<DEFNU>
% This syntax helper checks that model metadata fields are as expected.
% This is a MAPS check to ensure that users fill in model metadata for each
% model parsed in
%
% INPUTS:   
%   -> FileContents: structure containing all model info
%   -> FileKeywords: structure containing the keywords used in the file 
%   -> FileLineNumbers: structure containing the line numbers in the file
%   -> metadataFieldsCheckConfig: config for the metadata fields check
%
% OUTPUTS:  
%   -> none
%
% CALLS:
%   -> generate_MAPS_exception
%   -> generate_MAPS_exception_and_add_as_cause

%% UNPACK METADATA CONFIGURATION
% Unpack the file name of the metadata fields field used in the model file
% and the list of expected fields from the configuration.
metadataFieldsFieldName = metadataFieldsCheckConfig{1};
expectedMetadataFields = metadataFieldsCheckConfig{2};

%% CHECK FOR MISSING & UNKNOWN FIELDS
% If the model metadata fields are part of the model input, unpack them and
% compare them against the expected fields to find the logical index
% numbers of the missing & unknown fields.
if isfield(FileContents,metadataFieldsFieldName)
    metadataFields = FileContents.(metadataFieldsFieldName);
    missingFieldLogicals = find(cellfun(...
        @(x) ~any(strcmp(x,metadataFields)),expectedMetadataFields));
    unknownFieldLogicals = find(cellfun(...
        @(x) ~any(strcmp(x,expectedMetadataFields)),metadataFields));
end

%% CONSTRUCT AND THROW AN EXCEPTION
% If model metadata fields are part of the model input and there are
% missing or unknown fields, construct a master exception explaining the
% rule. Add causes to that exception for each of the two types of rule that
% were broken.
if isfield(FileContents,metadataFieldsFieldName)
   if any(missingFieldLogicals) || any(unknownFieldLogicals)
        masterErrId = ['MAPS:',mfilename,':BadMetadataFields']; 
        masterErrArgs = {FileKeywords.(metadataFieldsFieldName)};
        BadMetadataFieldsE = generate_MAPS_exception(...
            masterErrId,masterErrArgs);
        if any(missingFieldLogicals)
            errId = [masterErrId,':MissingFields'];
            errArgs{1} = expectedMetadataFields(missingFieldLogicals);
            BadMetadataFieldsE = ...
                generate_MAPS_exception_and_add_as_cause(...
                BadMetadataFieldsE,errId,errArgs);
        end
        if any(unknownFieldLogicals)
            metadataFieldLines = FileLineNumbers.(metadataFieldsFieldName);
            errId = [masterErrId,':UnknownFields'];
            errArgs{1} = metadataFields(unknownFieldLogicals)';
            errArgs{2} = ...
                num2str(metadataFieldLines(unknownFieldLogicals)');
            BadMetadataFieldsE = ...
                generate_MAPS_exception_and_add_as_cause(...
                BadMetadataFieldsE,errId,errArgs);          
        end
        throw(BadMetadataFieldsE);
   end   
end

end

%% FUNCTION TO CHECK THAT PARAMETER VALUES ARE VALID
function check_parameter_values_are_valid(...
    FileContents,FileKeywords,FileLineNumbers,parametersCheckConfig)        %#ok<DEFNU>
% This syntax helper checks that model parameters are proper numbers.
% This is a MAPS check to ensure that parameters parsed in (which are
% parsed) as strings can be converted properly to numeric values.
%
% INPUTS:   
%   -> FileContents: structure containing all model info
%   -> FileKeywords: structure containing the keywords used in the file 
%   -> FileLineNumbers: structure containing the line numbers in the file
%   -> parametersCheckConfig: config for the parameters check
%
% OUTPUTS:  
%   -> none

%% UNPACK CONFIGURATION INFO
% Unpack the parameter field name and the check type (either complete in
% which case all parameters must be provided or incomplete in which some
% can be left blank).
paramFieldName = parametersCheckConfig{1};
paramCheckType = parametersCheckConfig{2};

%% CONVERT CHECK TYPE TO BOOLEAN
% Convert the check type to a boolean to be used in the logical computation
% below.
if strcmp(paramCheckType,'complete')
    shouldParamsBeComplete = true;
else
    shouldParamsBeComplete = false;
end

%% CHECK PARAMETERS
% If parameters exist in the file contents, attempt a "regexp" match on
% each parameter. The expression allows for three alternative parameter 
% formats: '^(-)?((\d+)$' translates as "the entire string must be 
% comprised of optionally a minus sign, then a string of numbers";
% '^(-)?(\.)(\d+))$' translates as "the entire string must be comprised of
% optionally a minus sign, then a decimal point, then a string of numbers";
% and '^(-)?((\d+)(\.)(\d+))$'as  "the entire string must be comprised of
% optionally a minus sign, then a string of numbers, then a decimal point,
% then another string of numbers". For example, -0.3, 0.3, 01.3, .3 and 3.0
% are all valid numbers. If the checking config dictates that the list of 
% parameters can be incomplete, then a further check is made to remove any
% unmtached parameters in the regexp check that have been input as empty.
if isfield(FileContents,paramFieldName)
    paramStrs = FileContents.(paramFieldName);
    regexpNumExpr1 = '^(-)?(\d+)$';
    regexpNumExpr2 = '^(-)?(\.)(\d+)$';
    regexpNumExpr3 = '^(-)?(\d+)(\.)(\d+)$';
    regexpNumExpr = [regexpNumExpr1,'|',regexpNumExpr2,'|',regexpNumExpr3];
    matchedValidParams = regexp(paramStrs,regexpNumExpr,'match');
    badParamLogicals = cellfun(@isempty,matchedValidParams);
    if ~shouldParamsBeComplete
        badParamLogicals = (badParamLogicals&~cellfun(@isempty,paramStrs));
    end
    if any(badParamLogicals)
        paramKeyword = FileKeywords.(paramFieldName);
        paramLineNumbers = FileLineNumbers.(paramFieldName);
        masterErrId = ['MAPS:',mfilename,':BadParams'];
        masterErrArgs = {paramKeyword};
        errArgsList = [paramStrs ...
            convert_numeric_column_vector_to_string_equivalent(...
            paramLineNumbers)];
        generate_MAPS_exception_add_causes_and_throw(...
            masterErrId,errArgsList,badParamLogicals,masterErrArgs);
    end
end

end

%% HELPER FUNCTION TO UNPACK & COMBINE FILE CONTENTS
function [allContent,allContentLineNumbers,allContentKeywords] = ...
    unpack_and_combine_file_contents(...
    fieldNames,FileContents,FileLineNumbers,FileKeywords)
% This syntax check helper unpacks & combines specified model file fields.
% If passed the right number of inputs, it will also return the file line
% numbers and keywords describing all unpacked and combined content.
%
% INPUTS:   
%   -> fieldNames: field names of the file content to unpack & combine
%   -> FileContents: structure containing all model info
%   -> FileLineNumbers (optional): structure containing the line numbers
%   -> FileKeywords (optional): structure containing the keywords used 
%
% OUTPUTS:  
%   -> allContent: cell string array of all the fields' content
%   -> allContentLineNumbers: file line numbers corresponding to the
%      content
%   -> allContentKeywords: file keywords corresponding to the content
%
% CALLS:
%   -> none

%% SETUP OUTPUT
% Count the number of fields and prepare intermediate outputs (depending on
% the number of inputs passed in).
nFields = size(fieldNames,1);
unpackedContent = cell(1,nFields);
if nargin > 2
    unpackedContentLineNumbers = cell(1,nFields);
end
if nargin > 3
    unpackedContentKeywords = cell(1,nFields);
end

%% COLLECT ALL CONTENTS
% Run through each of the fields specified in the input and unpack them 
% into the cell array constructed above if they exist in the file contents. 
% Do the same for file line numbers and expand out the keyword related to 
% the field (depending on the number of inputs passed in).
for iField = 1:nFields
    iFieldName = fieldNames{iField};
    if isfield(FileContents,iFieldName)
        unpackedContent{iField} = FileContents.(iFieldName);
        if nargin > 2
            unpackedContentLineNumbers{iField} = ...
                FileLineNumbers.(iFieldName);
        end
        if nargin > 3
            unpackedContentKeywords{iField} = repmat(...
                {FileKeywords.(iFieldName)},size(unpackedContent{iField}));
        end
    end
end

%% CONCETENATE CONTENT INTO SINGLE VECTOR FOR OUTPUT
% Concatentate the content into a single vector for output using the MATLAB
% vertcat function.
allContent = vertcat(unpackedContent{:});
if isempty(allContent)
   allContent = {''}; 
end
if nargin > 2
    allContentLineNumbers = vertcat(unpackedContentLineNumbers{:});
end
if nargin > 3
    allContentKeywords = vertcat(unpackedContentKeywords{:});
end

end