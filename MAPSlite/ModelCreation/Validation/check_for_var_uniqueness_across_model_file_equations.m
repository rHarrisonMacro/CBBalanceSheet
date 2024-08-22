function EqSyntaxE = ...
    check_for_var_uniqueness_across_model_file_equations(...
    EqSyntaxE,FileContents,FileLineNumbers,eqsCheckConfig)
% This helper checks that there are not repeated variables across model
% file equations if such repetition is banned by the model file syntax
% rules.
%
% INPUTS:
%   -> EqSyntaxE: an exception to add causes to
%   -> FileContents: structure containing all the information in the model
%      file deconstructed into its constituent parts
%   -> FileKeywords: structure (with the same fields as the above) 
%      containing the keywords used in the model file
%   -> FileLineNumbers: structure (with the same fields as the above) 
%      containing the line numbers in the file from which the information
%      was taken 
%   -> eqsCheckConfig: cell array of configuration information for checks
%
% OUTPUTS:
%   -> EqSyntaxE: updated exception
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> is_row_cell_array
%   -> extract_LHS_expressions_from_equations
%   -> extract_RHS_expressions_from_equations
%   -> compute_equations_incidence_matrix
%   -> generate_MAPS_exception_and_add_as_cause
%   -> create_comma_separated_list
%
% DETAILS:
%   -> This helper is called by both the NLBL and LSS model file syntax
%      checkers.
%   -> It checks that a set of variables parsed from a model file are used
%      uniquely across a set of equations parsed from a model file.
%   -> It adds an exception if: (i) none of the variables appear on the 
%      specified side of the equation; (ii) more than one of the variables
%      appears on the specified side of the equation; (iii) a variable
%      appears in more than one equation (on the specified side).
%
% NOTES:
%   -> See <> for details of model file and add-on file syntax checking in
%      MAPS.
%
% This version: 15/06/2011
% Author(s): Alex Haberis & Matt Waldron

%% CHECK INPUTS
% Complete some basic checks on the input, including that the number of
% input arguments is correct.
if nargin < 4
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~strcmp(class(EqSyntaxE),'MException')
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(FileContents)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~isstruct(FileLineNumbers)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
elseif ~is_row_cell_array(eqsCheckConfig) || size(eqsCheckConfig,2)<3
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);    
end

%% UNPACK CONFIG INPUT
% Unpack the information needed for the checks: the name of the equations
% in the parsed file content, the name of the mnemonics and the side of 
% the equations to check.
typeOfEqsToCheck = eqsCheckConfig{1};
mnemsForUniqueVars = eqsCheckConfig{2};
eqSideForMnems = eqsCheckConfig{3};

%% UNPACK EQUATION STRINGS & LINE NUMBERS
% Unpack the equations and line numbers on which they appear. Exit the
% function if the equations are not part of the file contents
if ~isfield(FileContents,typeOfEqsToCheck)
    return
else
    eqStrs = FileContents.(typeOfEqsToCheck);
    eqLineNumbers = FileLineNumbers.(typeOfEqsToCheck);
end

%% UNPACK MNEMS FOR VARIABLES WHICH MUST APPEAR UNIQUELY
% Exit the function if the equations are not part of the file contents.
if ~isfield(FileContents,mnemsForUniqueVars)
    return
else
    varMnems = FileContents.(mnemsForUniqueVars);
end
    
%% DEFINE THE PART OF THE EQUATION IN WHICH VARIABLES MUST APPEAR UNIQUELY
% Either LHS, RHS or the whole equation
if strcmp(eqSideForMnems,'LHS')
    eqSidesToCheckStrs = extract_LHS_expressions_from_equations(eqStrs);
elseif strcmp(eqSideForMnems,'RHS')
    eqSidesToCheckStrs = extract_RHS_expressions_from_equations(eqStrs);
else
    eqSidesToCheckStrs = eqStrs;
end

%% COMPUTE INCIDENCE MATRIX
% Incidence matrix for variable appearances indexed by equations on rows
% and variables on columns.
eqVarIncMat = compute_equations_incidence_matrix(...
    eqSidesToCheckStrs,varMnems);

%% CHECK FOR EQUATIONS WITH NONE OF THE VARIABLES
% Find index numbers of equations with none of the variables. Construct
% exceptions for any such equations.
noVarEqInds = find(sum(eqVarIncMat,2)<1);
if ~isempty(noVarEqInds)
    errId = ['MAPS:',mfilename,':NoVarsInEq'];
    nBadEqs = size(noVarEqInds,1);
    for iEq = 1:nBadEqs
        iEqInd = noVarEqInds(iEq);
        errArgs = {num2str(eqLineNumbers(iEqInd))};
        EqSyntaxE = generate_MAPS_exception_and_add_as_cause(...
           EqSyntaxE,errId,errArgs);
    end
end

%% CHECK FOR EQUATIONS WITH TOO MANY OF THE VARIABLES
% Find index numbers of equations with more than one of the variables. 
% Construct exceptions for any such equations with the mnemonics of the
% variables found.
tooManyVarEqInds = find(sum(eqVarIncMat,2)>1);
if ~isempty(tooManyVarEqInds)
    errId = ['MAPS:',mfilename,':TooManyVarsInEq'];
    errArgs = cell(1,2);
    nBadEqs = size(tooManyVarEqInds,1);
    for iEq = 1:nBadEqs
        iEqInd = tooManyVarEqInds(iEq);
        errArgs{1} = num2str(eqLineNumbers(iEqInd));        
        errArgs{2} = create_comma_separated_list(...
            varMnems(eqVarIncMat(iEqInd,:))');
        EqSyntaxE = generate_MAPS_exception_and_add_as_cause(...
           EqSyntaxE,errId,errArgs);
    end
end

%% CHECK FOR REPETITIONS ACROSS EQUATIONS
% Find index numbers of variables which appear more than once in the 
% equations.  Construct exceptions for any such variables with the line
% numbers of the equations in which they appear.
badVarInds = find(sum(eqVarIncMat,1)>1)';
if ~isempty(badVarInds)
    errId = ['MAPS:',mfilename,':BadUseOfVarAcrossEqs'];
    errArgs = cell(1,2);
    nBadVars = size(badVarInds,1);
    for iVar = 1:nBadVars
        iVarInd = badVarInds(iVar);
        iVarEqLogicals = eqVarIncMat(:,iVarInd);        
        errArgs{1} = num2str(eqLineNumbers(iVarEqLogicals)');        
        errArgs{2} = varMnems{iVarInd};
        EqSyntaxE = generate_MAPS_exception_and_add_as_cause(...
           EqSyntaxE,errId,errArgs);
    end
end

end