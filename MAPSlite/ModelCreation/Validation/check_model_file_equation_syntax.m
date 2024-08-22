function EqSyntaxE = check_model_file_equation_syntax(EqSyntaxE,...
    FileContents,FileLineNumbers,eqCheckConfig,Model)
% This function validates the syntax of a set of model file equations.
% It checks that the equations are valid mathematical expressions in MATLAB
% and that their content conforms to the configuration information passed
% in.
%
% INPUTS:
%   -> EqSyntaxE: exception to add causes to
%   -> FileContents: structure containing all the information in the model
%      file deconstructed into its constituent parts
%   -> FileKeywords: structure (with the same fields as the above) 
%      containing the keywords used in the model file
%   -> FileLineNumbers: structure (with the same fields as the above) 
%      containing the line numbers in the file from which the information
%      was taken
%   -> eqCheckConfig: configuration information for the checks
%   -> Model: MAPS model object
%
% OUTPUTS:
%   -> EqSyntaxE: updated exception
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> unpack_model
%   -> remove_other_operators_from_equation (sub-function)
%   -> generate_MAPS_exception
%   -> check_equation_is_valid
%   -> check_equation_contents (sub-function)
%
% DETAILS:
%   -> This is a model file syntax check helper which checks the syntax of 
%      a generic set of modelf file equations according to a configuration. 
%   -> If it finds any syntax errors, it adds the exceptions as causes to
%      the input exception passed in.
%   -> The syntax checks are split into two parts: those that relate to
%      the validity of terms and expressions in MATLAB (i.e. are the
%      expressions valid MATLAB expressions and are all terms valid MATLAB
%      variable names) and those that relate to the specifics of equation
%      content given the rest of the information in the model file (and
%      possibly in the model itself).
%   -> The first check is largely configuration free, while the second
%      checks that the terms and time subscripts used in each of the 
%      equations conforms to the information in the configuration and in 
%      the file content. 
%   -> For example, there are two model equations in the NKPC LSS model.
%      This function would first check that those equations are valid (eg
%      x{t} = rho*x{t-1} + sx*ux{t} is valid but x{t} = rho*/x{t-1} +
%      sx*ux{t} is not). It would then check that the equations contain
%      only model variables, shocks, and parameters and that the time
%      subscripts associated with those are valid.
%   -> This function also deals with checking of equations in a file that 
%      should be associated with an existing model (through the final model
%      input) such as those in add-on files. In such cases, the
%      congfiguration information also describes whether the equation
%      content is expected to be found in the file contents or the model
%      and, if in the model, the model characteristics under which it is
%      expected to be found.
%
% NOTES:
%   -> See <> for more details of MAPS model creation.
%   -> See also <> for a description of MAPS model file syntax checking
%      rules.
%
% This version: 28/04/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of inputs is as expected.
if nargin < 4
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});  
elseif ~strcmp(class(EqSyntaxE),'MException');
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(FileContents)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);   
elseif ~isstruct(FileLineNumbers)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);  
elseif ~is_row_cell_array(eqCheckConfig) || size(eqCheckConfig,2)<3
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);   
elseif nargin > 4 
    if ~isstruct(Model)
        errId = ['MAPS:',mfilename,':BadInput5'];
        generate_and_throw_MAPS_exception(errId);
    end
end

%% UNPACK CONFIG INPUT
% The configuration input is formatted such that it is a 1*3/4 cell
% array with the 1st element containing the string field name of the
% equations to check, the 2nd element revealing the type of equation
% ('explicit' or 'implicit') and the 3rd containing information about the
% permitted content (field names of the variables / parameters, time
% subscripts that can be associated with those and information aboput
% whether they should only appear on the left-hand-side, right-hand-side or
% non-specific). Finally, the optional 4th element describes additional 
% operators (above and beyond those that are listed in the valid MAPS 
% operators config that may appear in the equation. If a 5th model input was 
% passed in unpack the additional configuration information which describes
% whether the expected equation content exists in the file content or the
% model and the condition under which it lives in the model (eg if 
% modelHasSteadyStateEqs).
fieldNameOfEqsTocheck = eqCheckConfig{1};
typeOfEqsToCheck = eqCheckConfig{2};
permittedEqContentConfig = eqCheckConfig{3};
fieldNamesOfMnems = permittedEqContentConfig(:,1);
timeSubsForMnems = permittedEqContentConfig(:,2);
eqSideForMnems = permittedEqContentConfig(:,3);
if nargin > 4
    locationOfMnems = permittedEqContentConfig(:,4);
    conditionForMnems = permittedEqContentConfig(:,5);
end
if size(eqCheckConfig,2) > 3
    confgContainsAdditionalOps = true;
    otherPermittedOps = eqCheckConfig{4};
else
    confgContainsAdditionalOps = false;
end

%% UNPACK EQUATION STRINGS & LINE NUMBERS
% Unpack the equations and line numbers on which they appear.
eqStrs = FileContents.(fieldNameOfEqsTocheck);
eqLineNumbers = FileLineNumbers.(fieldNameOfEqsTocheck);

%% CREATE LIST OF PERMITTED TERMS WITH THEIR TIME SUBSCRIPTS
% Create a list of the mnemonics that are allowed to appear in the
% equations, the time subscripts for those mnemonics and the side of the
% equation on which they should appear. Note that this code makes no
% assumption about which fields should exist in the parsed model file 
% contents.
nMnemTypes = size(fieldNamesOfMnems,1);
mnemsCell = cell(1,nMnemTypes);
mnemsTimeSubsCell = cell(1,nMnemTypes);
mnemsEqSidesCell = cell(1,nMnemTypes);
for iType = 1:nMnemTypes
    if nargin < 5 || strcmp(locationOfMnems{iType},'file')
        if isfield(FileContents,fieldNamesOfMnems{iType})
            mnemsCell{iType} = FileContents.(fieldNamesOfMnems{iType});
        end
    elseif strcmp(locationOfMnems{iType},'model')
       if isempty(conditionForMnems{iType}) || ...
               unpack_model(Model,conditionForMnems(iType))
           mnemsCell{iType} = unpack_model(Model,fieldNamesOfMnems(iType));           
       end        
    end 
    niTypeMnems = size(mnemsCell{iType},1);
    mnemsTimeSubsCell{iType} = repmat(...
        timeSubsForMnems(iType),[niTypeMnems,1]);
    mnemsEqSidesCell{iType} = repmat(...
        eqSideForMnems(iType),[niTypeMnems,1]);
end
permittedTerms = vertcat(mnemsCell{:});
permittedTimeSubs = vertcat(mnemsTimeSubsCell{:});
permittedEqSides = vertcat(mnemsEqSidesCell{:});

%% AMMEND LIST FOR OTHER PERMITTED OPERATORS (IF APPLICABLE)
% If any other additional operators were passed in, add them to the list of
% valid terms (with empty time subscripts and non-specific equation sides).
if confgContainsAdditionalOps
    permittedTerms = [permittedTerms;otherPermittedOps];
    nOtherPermittedOps = size(otherPermittedOps,1);
    permittedTimeSubs = [permittedTimeSubs;...
        repmat({''},[nOtherPermittedOps,1])];
    permittedEqSides = [permittedEqSides;...
        repmat({''},[nOtherPermittedOps,1])];   
end

%% AMEND EQUATIONS FOR VALIDITY CHECK
% If the equations may contain additional operators, create a new set of
% equations with those additional operators removed and replaced with the
% 'log' operator. This is a workaround to problems associatd with the 
% "diff" operator used in LSS model data transformations which is (rightly)
% not recognised as a valid operator.
if confgContainsAdditionalOps
    eqStrsForValidityCheck = remove_other_operators_from_equation(...
        eqStrs,otherPermittedOps);
else
    eqStrsForValidityCheck = eqStrs;
end

%% SET LOGICAL EQUATION TYPE FIELD
% Use the explicit/implicit equation type field to set a logical variable
% for use in the equation validity checking function.
if strcmp(typeOfEqsToCheck,'explicit')
    eqTypeIsExplicit = true;
elseif strcmp(typeOfEqsToCheck,'implicit')
    eqTypeIsExplicit = false;
end

%% CHECK EQUATION IS A VALID EQUATION AND HAS THE RIGHT CONTENT
% Run through the equations one by one, checking their validity and then
% their content. If any exceptions are found, they are added as cause to
% the (equation-specific) exception created, which is then added as cause
% to the master exception passed in as input.
nEqStrs = size(eqStrs,1);
for iEq = 1:nEqStrs
    errId = ['MAPS:',mfilename,':BadEquation'];
    errArgs = {num2str(eqLineNumbers(iEq)),eqStrs{iEq}};
    EqSyntaxEi = generate_MAPS_exception(errId,errArgs);
    EqSyntaxEi = check_equation_is_valid(...
        eqStrsForValidityCheck{iEq},eqTypeIsExplicit,EqSyntaxEi);
    EqSyntaxEi = check_equation_contents(EqSyntaxEi,eqStrs{iEq},...
        permittedTerms,permittedTimeSubs,permittedEqSides);
    if ~isempty(EqSyntaxEi.cause)
        EqSyntaxE = addCause(EqSyntaxE,EqSyntaxEi);
    end    
end

end

%% FUNCTION TO CHECK EQUATION CONTENTS
function EqSyntaxEi = check_equation_contents(EqSyntaxEi,...
    eqStr,permittedTerms,permittedTimeSubs,permittedEqSides)
% This function validates the content of a model file equation.
% It checks that all terms in the equation are recognised, have the correct
% time subscripts and appear on the expected side of the equation.
%
% INPUTS:
%   -> EqSyntaxEi: equation exception to add causes to
%   -> eqStr: equation string
%   -> permittedTerms: list of permitted terms (other than numbers etc)
%   -> permittedTimeSubs: list of time subscripts permitted for each term
%   -> permittedEqSides: list of sides of the equation permitted for each 
%      term 
%
% OUTPUTS:
%   -> EqSyntaxEi: updated equation exception
%
% CALLS:
%   -> split_equation
%   -> generate_MAPS_exception_and_add_as_cause

%% SPLIT EQUATION INTO TERMS AND TIME SUBSCRIPTS
% Use the symbolic MAPS equation splitter to split the equation into terms,
% recognised delimiters (eg '+', '=' etc) and time subscripts associated
% with each of the terms.
[eqStrTerms,eqStrDelims,eqStrTimeSubs] = split_equation(eqStr);

%% FIND POSITION OF THE EQUAL SYMBOL
% Find the position of the equal sign in the equation string delimiters. If
% it does not exist, set the position equal zero (meaning all terms will be
% assumed on the LHS).
if any(strcmp('=',eqStrDelims))
    eqEqualSignPos = find(strcmp('=',eqStrDelims),1);
else
    eqEqualSignPos = 0;
end

%% FIND RECOGNISED & EMPTY TERMS
% Compute the number of terms found. For each, compute a logical vector
% which describes whether the term is empty, is recognised in the list of
% permitted terms and appears before the equal sign.
nTerms = size(eqStrTerms,2);
emptyTermLogicals = false(1,nTerms);
recognisedTermLogicals = false(1,nTerms);
beforeEqualSignTermLogicals = false(1,nTerms);
for iTerm = 1:nTerms
    emptyTermLogicals(iTerm) = isempty(eqStrTerms{iTerm});
    recognisedTermLogicals(iTerm) = any(strcmp(...
        eqStrTerms{iTerm},permittedTerms));
    beforeEqualSignTermLogicals(iTerm) = (iTerm<=eqEqualSignPos);
end

%% GATHER INFORMATION ABOUT RECOGNISED TERMS IN EQUATION STRING
% Use the logicals to compute the list of valid terms found, the time
% subscripts associated with them and whether or not they appear before the
% equal sign.
validTerms = eqStrTerms(recognisedTermLogicals)';
validTermsTimeSubs = eqStrTimeSubs(recognisedTermLogicals)';
validTermsBeforeEqualSign = beforeEqualSignTermLogicals(...
    recognisedTermLogicals)';

%% GATHER DETAILS ABOUT RECOGNISED TERMS
% Find the index numbers of the valid terms in the permitted terms and find
% the time subscripts and equation sides on which they are allowed to
% appear.
nValidTerms = size(validTerms,1);
indsOfValidTermsInList = zeros(nValidTerms,1);
for iTerm = 1:nValidTerms
    indsOfValidTermsInList(iTerm) = find(strcmp(...
        validTerms{iTerm},permittedTerms),1);
end
validTermsPermittedTimeSubs = permittedTimeSubs(indsOfValidTermsInList);
validTermsEqSides = permittedEqSides(indsOfValidTermsInList);

%% FIND UNRECOGNISED TERMS
% Find any terms that were both unrecognised & non-empty. If any such terms
% exist add an exception as cause with details of which tersm were
% unrecognised.
unrecognisedTermLogicals = (~recognisedTermLogicals&~emptyTermLogicals);
if any(unrecognisedTermLogicals)
    unrecognisedTerms = unique(eqStrTerms(unrecognisedTermLogicals));
    errId = ['MAPS:',mfilename,':UnrecognisedTerms'];
    errArgs = unrecognisedTerms;
    EqSyntaxEi = generate_MAPS_exception_and_add_as_cause(...
        EqSyntaxEi,errId,errArgs);
end

%% FIND TERMS WITH INVALID TIME SUBSCRIPTS
% For each of the valid terms, check that the time subscripts associated
% with them are valid. If not add exceptions as causes.
for iTerm = 1:nValidTerms
    iTermTimeSubs = validTermsTimeSubs{iTerm};
    iTermPermittedTimeSubs = validTermsPermittedTimeSubs{iTerm};
    if ~isempty(iTermTimeSubs) || ~isempty(iTermPermittedTimeSubs)        
        if isempty(regexp(iTermTimeSubs,iTermPermittedTimeSubs,'match'))
            errId = ['MAPS:',mfilename,':BadTimeSubscript'];
            errArgs = {validTerms{iTerm} iTermTimeSubs};
            EqSyntaxEi = generate_MAPS_exception_and_add_as_cause(...
                EqSyntaxEi,errId,errArgs);
        end
    end    
end

%% FIND TERMS ON WRONG SIDE OF EQUATION
% For each of the valid terms, check that they appear on the correct side
% of the equation (if one is specified). If not, add exceptions as cause.
for iTerm = 1:nValidTerms
    if (strcmp(validTermsEqSides{iTerm},'LHS')&&...
            ~validTermsBeforeEqualSign(iTerm)) || ...
            (strcmp(validTermsEqSides{iTerm},'RHS')&&...
            validTermsBeforeEqualSign(iTerm))
        errId = ['MAPS:',mfilename,':WrongEquationSide'];
        errArgs = {validTerms{iTerm} validTermsEqSides{iTerm}};
        EqSyntaxEi = generate_MAPS_exception_and_add_as_cause(...
            EqSyntaxEi,errId,errArgs);
    end
end

end

%% FUNCTION TO REMOVE ADDITIONAL PERMITTED (NON-MATLAB) OPERATORS
function eqStrsForValidityCheck = remove_other_operators_from_equation(...
    eqStrs,otherPermittedOps)
% This function removes addition / unrecognised operators from equations.
% It is built as a workaround to the "diff" operator used in LSS data
% transformations, which causes a problem in the equation validity check
% (because MATLAB righly recognises "diff" as a time series operator).
%
% INPUTS:
%   -> eqStrs: set of equation strings
%   -> otherPermittedOps: cell string array list of "special" operators 
%
% OUTPUTS:
%   -> eqStrsForValidityCheck: amended set of equation strings
%
% CALLS:
%   -> split_equation
%   -> reconstruct_equation    

%% SETUP OUTPUT
% Compute the number of equations input and setup the output.
nEqs = size(eqStrs,1);
eqStrsForValidityCheck = cell(nEqs,1);

%% REPLACE ALL SEPICAL OPERATORS WITH LOG
% Replace all instances of special operators found in all the equations 
% with the 'log' operator which is known to be valid. 
nOtherPermittedOps = size(otherPermittedOps,1);
for iEq = 1:nEqs
    [ieqStrTerms,ieqStrDelims,ieqStrTimeSubs] = split_equation(...
        eqStrs{iEq});
    for iOp = 1:nOtherPermittedOps
        ieqStrTerms(strcmp(otherPermittedOps{iOp},ieqStrTerms)) = {'log'};
    end
    eqStrsForValidityCheck{iEq} = reconstruct_equation(...
        ieqStrTerms,ieqStrDelims,ieqStrTimeSubs);
end

end