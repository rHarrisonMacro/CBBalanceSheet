function ModelwTgtRule = ...
    create_LSS_model_with_discretion_targeting_rule(Model,ODinfo,...
    tgtRuleNames,useModifierToCreateModel,tgtRuleModifierFileName,...
    paramsToDiscardMnems,OptPolSolOptions)
% Creates a MAPS LSS model with an optimal discretion targeting rule.
% The targeting rule is computed using a MAPS implementation of the Dennis
% algorithm. The model is then solved and created using the standard MAPS
% toolkit which uses AIM.
%
% INPUTS:   
%   -> Model: MAPS LSS model
%   -> ODinfo: structure of info required to solve the OD problem
%       - policyEqNames: single string or column cell string array
%       - policyShockMnems: single string or column cell string array (use
%         emmpty string ('') if there are no policy shocks to remove)
%       - instrumentMnems: single string or column cell string array
%       - instrumentWeights: column vector of loss weights for instruments
%       - objVarMnems: single string or column cell string array
%       - objVarWeights: column vector of loss weights for objective vars
%       - beta: discount factor
%   -> tgtRuleNames: string or column cell string array of names for the
%      targeting rules
%   -> useModifierToCreateModel: false/true
%   -> tgtRuleModifierFileName: '.maps' extended file name for the
%      discretion targeting rule(s) modifer file (use an empty string ('')
%      if useModifierToCreateModel=false)
%   -> paramsToDiscardMnems: names of any parameters to remove
%      from the model (use an empty string ('') if there are none to 
%      remove)
%   -> OptPolSolOptions (optional): structure
%       - algorithmIsDennis (optional): true/false
%       - Initialisation (optional): Initialisation choice for B/V
%           - source (optional): 'setToZeros', 'useInputModel', 
%             'setToUserInput'
%           - Binit: nx*nx matrix for `setToUserInput' option
%           - VxtildextildeInit: nxtilde*nxtilde matrix for 
%             `setToUserInput' option if algorithmIsDennis = false
%       - Algorithm (optional): structure
%           - tol (optional): tolerance for convergence
%           - maxIter (optional): maximum number of iterations
%
% OUTPUTS:  
%   -> ModelwTgtRule: model with targeting rule(s) in place of existing
%      policy rules
%
% DETAILS: 
%   -> This function creates a MAPS model with a discretion targeting rule
%      (or multiple rules) in it.
%   -> The discretion targeting rule is computed using a MAPS
%      implementation of the Dennis algorithm (or VFI depending on the 
%      Options). See ANALYTICAL 7702486 for a full derivation.
%   -> In addition to the computation of the targeting rule, this function
%      also amends the model by removing the existing policy shocks and any
%      paramaters that formed part of the now-redundant policy rules.
%   -> It is possible that the model will exhibit indeterminacy in the
%      sense that the targeting rule(s) may not pin down a unique RE
%      solution.  In that case, the modification process will fail and the
%      code below will construct the model "by hand".
%   -> The model that emerges from this function is not compatible with
%      the use of anticipated shocks or with imposing bound constraints on
%      the instruments via the Holden and Paetz methodlogy: reflecting this
%      invalidity, the forward loadings matrix is set to zero.
%
% NOTES:
%   -> The precision on the loadings in the targeting rule is reduced in
%      the code below - the intention of this is to eliminate loadings that
%      should in fact be 0 and to make the resulting targeting rule
%      equation string more readable.
%   -> If the modification process is successful, then the targeting rules
%      are appended to the end of the model.  To the extent that it makes 
%      any sense to associate the names of the targeting rules with 
%      particular instruments, the ordering of the targeting rule names 
%      input to this function should be made consistent with the ordering 
%      of the instrument mnemonics in the Odinfo structure input.
%   -> Error handling is kept to a minimum with some of it "delegated" to
%      functions called below.
%
% This version: 23/08/2018
% Author(s): Matt Waldron (based on initial code by Alex Haberis)

%% CHECK INPUTS
if nargin < 6
    error('This function requires 6 or 7 inputs')
end

%% HANDLE OPTIONAL ALGORITHM OPTIONS INPUT
if nargin < 7
    OptPolSolOptions = struct;
end

%% SOLVE MODEL FOR OPTIMAL DISCRETIONARY POLICY
[BOD,PHIOD,HBOD,HCOD,HFOD,PSIOD,~,ODinfo] = ...
    solve_LSS_model_under_optimal_discretion(...
    Model,ODinfo,OptPolSolOptions);

%% CHECK CONSISTENCY OF TARGETING RULE NAME DIMENSION WITH POLICY EQ DIM
policyEqNames = ODinfo.policyEqNames;
if ~is_string_or_column_cell_string_array(tgtRuleNames) || ...
        size(tgtRuleNames,1) ~= size(policyEqNames,1)
    error(['3rd input must be a string or column cell string array ',...
        'of names for the targeting rule equation(s), and must match ',...
        'the number of policy equation names, which they are ',...
        'replacing (and which form part of the ODinfo structure)']);
end

%% CHECK PARAMS TO DISCARD INPUT
if ~is_string_or_column_cell_string_array(paramsToDiscardMnems)
    error('6th input must be a string or column cell string array');
end

%% HOMEGENISE TARGETING/POLICY EQ NAMES & PARAMS TO DISCARD MNEMS TO CELLS
% This is just so that the code below doesn't have to take care of the 
% alternative input formats.
nTgtRules = size(tgtRuleNames,1);
if nTgtRules==1 
    if ischar(tgtRuleNames)
        tgtRuleNames = {tgtRuleNames};
    end
    if ischar(policyEqNames)
        policyEqNames = {policyEqNames};
    end
end
if ischar(paramsToDiscardMnems)
    if isempty(paramsToDiscardMnems)
        paramsToDiscardMnems = {};
    else
        paramsToDiscardMnems = {paramsToDiscardMnems};
    end
end
ODinfo.policyEqNames = policyEqNames;

%% CREATE MODEL USING MODIFIER (IF RELEVANT)
if useModifierToCreateModel
    try
        ModelwTgtRule = create_OD_model_using_modifier(...
            Model,ODinfo,tgtRuleNames,tgtRuleModifierFileName,...
            paramsToDiscardMnems,HCOD);
    catch ModelCreationE
        warning(['Model modifier function failed with the following ',...
            'message ''',ModelCreationE.message, ''' ... constructing ',...
            'the OD model ''by hand'' instead']);
        ModelwTgtRule = create_OD_model_by_hand(...
            Model,ODinfo,tgtRuleNames,paramsToDiscardMnems,BOD,PHIOD,...
            HBOD,HCOD,HFOD,PSIOD);
    end
end

%% CREATE MODEL BY HAND (IF RELEVANT)
if ~useModifierToCreateModel
    ModelwTgtRule = create_OD_model_by_hand(Model,ODinfo,tgtRuleNames,...
        paramsToDiscardMnems,BOD,PHIOD,HBOD,HCOD,HFOD,PSIOD);
end

%% OVERWRITE F MATRIX WITH ZEROS
nx = size(BOD,1);
F = zeros(nx,nx);
ModelwTgtRule = pack_model(ModelwTgtRule,'F',{F});

%% INSERT CUSTOM SOLVER
% Note that the "fixed" inputs in the function handle are used to recompute
% the OD solution given the new set of parameters have been adapted to
% reflect the amendments to the model given above.  In particular, the
% policy shocks have been removed (and so the policyShockMnems
% input is set to empty) and the policy rules have become targeting
% equations.
ODinfoForResolve = ODinfo;
ODinfoForResolve.policyShockMnems = '';
ODinfoForResolve.policyEqNames = tgtRuleNames;
ModelwTgtRule.Numerics.customSolver = ...
    @(Model,paramVec) resolve_LSS_model_with_discretion_targeting_rule(...
    Model,paramVec,ODinfoForResolve,OptPolSolOptions);

end

%% FUNCTION TO CREATE OD MODEL USING MODIFIER
function ModelwTgtRule = create_OD_model_using_modifier(Model,ODinfo,...
    tgtRuleNames,tgtRuleModifierFileName,paramsToDiscardMnems,HCOD)
% This function uses the model modifier to create a model with a tgt rules.
% The model is then solved using AIM.  It may not reproduce the OD solution
% exactly and/or maybe indeterminate.
%
% INPUTS:   
%   -> Model: MAPS LSS model
%   -> ODinfo: structure of info required to solve the OD problem
%       - policyEqNames: single string or column cell string array
%       - policyShockMnems: single string or column cell string array (use
%         emmpty string ('') if there are no policy shocks to remove)
%       - instrumentMnems: single string or column cell string array
%       - instrumentWeights: column vector of loss weights for instruments
%       - objVarMnems: single string or column cell string array
%       - objVarWeights: column vector of loss weights for objective vars
%       - beta: discount factor
%       - policyEqLogicals: vector of logicals that locate policy eqs
%       - policyEqInds: vector of indices that locate policy eqs
%       - policyShockLogicals: vector of logicals locating policy shocks
%       - policyShockInds: vector of indices locating policy shocks
%       - instrumentLogicals: vector of logicals locating the instruments
%       - instrumentInds: vector of indices that locate the instruments
%       - Q: instruments weights matrix
%       - W: target weights matrix
%   -> tgtRuleNames: column cell string array of names for the tgt rules
%   -> tgtRuleModifierFileName: '.maps' extended file name for the
%      discretion targeting rule(s) modifer file (use an empty string ('')
%      if useModifierToCreateModel=false)
%   -> paramsToDiscardMnems: column cell string array for the parameters to
%      discard
%   -> HCOD: loadings on contemporaneous variables in structural eqs
%
% OUTPUTS:  
%   -> ModelwTgtRule: model with targeting rule(s) in place of existing
%      policy rules

%% CONSTRUCT POLICY RULES REMOVAL INSTRUCTION CELL
policyRuleRemovalsCell = strcat('Remove:',ODinfo.policyEqNames);

%% CONSTRUCT TARGETING RULE ADDITION INSTRUCTION CELL
nDigitsAfterDecimal = 12;
nSignificantFigs = 12;
tgtRuleEqStrs = construct_targeting_rule_equation_strings(Model,ODinfo,...
    HCOD,nDigitsAfterDecimal,nSignificantFigs);
tgtRuleAdditionsCell = strcat('Add:',tgtRuleNames,':',tgtRuleEqStrs);

%% CONSTRUCT INSTRUCTION TO DELETE REDUNDANT POLICY SHOCKS
policyShockMnems = ODinfo.policyShockMnems;
nShocksToDiscard = size(policyShockMnems,1);
if nShocksToDiscard > 0
    shocksRemovalCell = strcat('Remove:',policyShockMnems);
end
    
%% CONSTRUCT INSTRUCTION TO DELETE REDUNDANT POLICY RULE PARAMATERS
nParamsToDiscard = size(paramsToDiscardMnems,1);
if nParamsToDiscard > 0
    paramsRemovalCell = strcat('Remove:',paramsToDiscardMnems);
end

%% CREATE MODIFIER FOR DISCRETION TARGETING RULE
contentToWrite = [{...
    'METADATA'
    'Name: createODtgtRuleModel'
    ['Description: Modify monetary policy rule to implement targeting ',...
    'rule for discretion']
    'Author: autogenerated'
    ''
    'MODEL EQUATIONS'
    }
    tgtRuleAdditionsCell;
    policyRuleRemovalsCell];
if nShocksToDiscard > 0
    contentToWrite = [contentToWrite;{'';'SHOCKS'};shocksRemovalCell];
end
if nParamsToDiscard > 0
    contentToWrite = [contentToWrite;{'';'PARAMETERS'};paramsRemovalCell];
end

%% WRITE MODIFICATION INFO TO FILE
write_to_text_file(tgtRuleModifierFileName,contentToWrite);

%% MODIFY MODEL
ModelwTgtRule = modify_model_from_struct(Model,tgtRuleModifierFileName);

end

%% FUNCTION TO CREATE OD MODEL BY HAND
function ModelwTgtRule = create_OD_model_by_hand(Model,ODinfo,...
    tgtRuleNames,paramsToDiscardMnems,BOD,PHIOD,HBOD,HCOD,HFOD,PSIOD)
% This function reproduces a modification using the modifier.
% NB. The only difference is that this function retains the ordering of the
% original policy rules.
%
% INPUTS:   
%   -> Model: MAPS LSS model
%   -> ODinfo: structure of info required to solve the OD problem
%       - policyEqNames: single string or column cell string array
%       - policyShockMnems: single string or column cell string array (use
%         emmpty string ('') if there are no policy shocks to remove)
%       - instrumentMnems: single string or column cell string array
%       - instrumentWeights: column vector of loss weights for instruments
%       - objVarMnems: single string or column cell string array
%       - objVarWeights: column vector of loss weights for objective vars
%       - beta: discount factor
%       - policyEqLogicals: vector of logicals that locate policy eqs
%       - policyEqInds: vector of indices that locate policy eqs
%       - policyShockLogicals: vector of logicals locating policy shocks
%       - policyShockInds: vector of indices locating policy shocks
%       - instrumentLogicals: vector of logicals locating the instruments
%       - instrumentInds: vector of indices that locate the instruments
%       - Q: instruments weights matrix
%       - W: target weights matrix
%   -> tgtRuleNames: column cell string array of names for the tgt rules
%   -> paramsToDiscardMnems: column cell string array for the parameters to
%      discard
%   -> BOD: state transition for model solved under optimal discretion
%   -> PHIOD: loadings on shocks in optimal discretion solution
%   -> HBOD: loadings on lagged variables in structural eqs
%   -> HCOD: loadings on contemporaneous variables in structural eqs
%   -> HFOD: loadings on one-period ahead variables in structural eqs
%   -> PSIOD: loadings on shocks in structural eqs
%
% OUTPUTS:  
%   -> ModelwTgtRule: model with targeting rule(s) in place of existing
%      policy rules

%% UNPACK EQUATION NAMES, PARAMETERS, PARAMETER MNEMS AND SHOCK MNEMS
[xEqODstrs,xEqODnames,theta,thetaMnems,thetaNames,zMnems,zNames] = ...
    unpack_model(Model,{'xEqStrs','xEqNames','theta','thetaMnems',...
    'thetaNames','zMnems','zNames'});

%% CONSTRUCT TARGETING RULE STRING EQUATIONS
nDigitsAfterDecimal = 12;
nSignificantFigs = 12;
tgtRuleEqStrs = construct_targeting_rule_equation_strings(Model,ODinfo,...
    HCOD,nDigitsAfterDecimal,nSignificantFigs);

%% OVERWRITE POLICY RULE EQUATION STRINGS & NAMES WITH TARGETING RULE NAMES
xEqODnames(ODinfo.policyEqInds) = tgtRuleNames;
xEqODstrs(ODinfo.policyEqInds) = tgtRuleEqStrs;

%% TAKE CARE OF PARAMETERS TO REMOVE
thetaODlogicals = ~ismember(thetaMnems,paramsToDiscardMnems);
thetaOD = theta(thetaODlogicals);
thetaODmnems = thetaMnems(thetaODlogicals);
thetaODnames = thetaNames(thetaODlogicals);

%% REMOVE POLICY SHOCKS            
PSIODnoPolicyShocks = PSIOD(:,~ODinfo.policyShockLogicals);
PHIODnoPolicyShocks = PHIOD(:,~ODinfo.policyShockLogicals);
zODnoPolicyShocksMnems = zMnems(~ODinfo.policyShockLogicals);
zODnoPolicyShocksNames = zNames(~ODinfo.policyShockLogicals);

%% COMPUTE VARIANCE COVARIANCE MATRIX
POD = compute_model_variable_covariance_matrix(PHIOD,BOD);

%% UPDATE MODEL
ModelwTgtRule = Model;
ModelwTgtRule = pack_model(ModelwTgtRule,{'B','PHI','P','HB','HC','HF',...
    'PSI','xEqStrs','xEqNames','theta','thetaMnems','thetaNames',...
    'zMnems','zNames'},{BOD,PHIODnoPolicyShocks,POD,HBOD,HCOD,HFOD,...
    PSIODnoPolicyShocks,xEqODstrs,xEqODnames,thetaOD,thetaODmnems,...
    thetaODnames,zODnoPolicyShocksMnems,zODnoPolicyShocksNames});

%% UPDATE MODEL SYMBOLICS
ModelwTgtRule = create_LSS_model_symbolics(ModelwTgtRule);

end

%% FUNCTION TO CREATE TARGETING RULE STRINGS
function tgtRuleEqStrs = construct_targeting_rule_equation_strings(...
    Model,ODinfo,HCOD,nDigitsAfterDecimal,nSignificantFigs)
% This function creates the targeting rules as string equations.
%
% INPUTS:   
%   -> Model: MAPS LSS model
%   -> ODinfo: structure of info required to solve the OD problem
%       - policyEqNames: single string or column cell string array
%       - policyShockMnems: single string or column cell string array (use
%         emmpty string ('') if there are no policy shocks to remove)
%       - instrumentMnems: single string or column cell string array
%       - instrumentWeights: column vector of loss weights for instruments
%       - objVarMnems: single string or column cell string array
%       - objVarWeights: column vector of loss weights for objective vars
%       - beta: discount factor
%       - policyEqLogicals: vector of logicals that locate policy eqs
%       - policyEqInds: vector of indices that locate policy eqs
%       - policyShockLogicals: vector of logicals locating policy shocks
%       - policyShockInds: vector of indices locating policy shocks
%       - instrumentLogicals: vector of logicals locating the instruments
%       - instrumentInds: vector of indices that locate the instruments
%       - Q: instruments weights matrix
%       - W: target weights matrix
%   -> HCOD: loadings on contemporaneous variables in structural eqs
%   -> nDigitsAfterDecimal: non-negative scalar
%   -> nSignificantFigs: non-negative scalar
%
% OUTPUTS:  
%   -> tgtRuleEqStrs: nRules*1 cell array of equation strings

%% UNPACK MODEL VARIABLE MNEMONICS
xMnems  = unpack_model(Model,'xMnems');

%% EXTRACT TARGETING RULE COEFFICIENTS
policyEqInds = ODinfo.policyEqInds;
tgtRuleCoeffs = HCOD(policyEqInds,:);

%% CONSTRUCT STRING EQUATION
nTgtRules = size(tgtRuleCoeffs,1);
tgtRuleEqStrs = cell(nTgtRules,1);
for iRule = 1:nTgtRules
    [iTgtRuleSignifCoeffStrs,iTgtRuleCoeffsNonZeroInd] = ...
        convert_numeric_coefficients_to_strings(tgtRuleCoeffs(iRule,:)',...
        nDigitsAfterDecimal,nSignificantFigs);        
    iTgtRuleSignifCoeffMnems = xMnems(iTgtRuleCoeffsNonZeroInd,1);
    nVarsToInclude = size(iTgtRuleSignifCoeffMnems,1);
    iTgtRuleStr = '';
    for iVarToInclude = 1:nVarsToInclude
        eqStrTerm = [iTgtRuleSignifCoeffStrs{iVarToInclude}...
            ,'*',iTgtRuleSignifCoeffMnems{iVarToInclude},'{t}'];
        if iVarToInclude < nVarsToInclude
            iTgtRuleStr = [iTgtRuleStr,eqStrTerm,'+'];                      %#ok<AGROW>
        else
            iTgtRuleStr = [iTgtRuleStr,eqStrTerm];                          %#ok<AGROW>
            iTgtRuleStr = [iTgtRuleStr,' = 0'];                             %#ok<AGROW>
            iTgtRuleStr = strrep(iTgtRuleStr,'+-','-');
        end
    end
    tgtRuleEqStrs{iRule} = iTgtRuleStr;
end

end

%% FUNCTION TO CONVERT NUMERIC COEFFICIENTS TO STRINGS
function [iTgtRuleSignifCoeffStrs,iTgtRuleCoeffsNonZeroInd] = ...
        convert_numeric_coefficients_to_strings(iTgtRuleCoeffs,...
        nDigitsAfterDecimal,nSignificantFigs)
% Converts the numeric vector of coefficients to strings.
%
% INPUTS:   
%   -> iTgtRuleCoeffs: nx*1 vector of numeric coefficients
%   -> nDigitsAfterDecimal: number of digits to allow after decimal place
%   -> nSignificantFigs: number of significant figures (for coefficients>1)
%
% OUTPUTS:  
%   -> iTgtRuleSignifCoeffStrs: nSignifCoeffs*1 column cell string array 
%   -> iTgtRuleCoeffsNonZeroInd: index numbers for the associated variables

%% FIND ALL COEFFICIENTS BIGGER THAN THE MINIMUM
iTgtRuleCoeffsNonZeroInd = (abs(iTgtRuleCoeffs)>...
    (1/(10^nDigitsAfterDecimal)));
iTgtRuleCoeffsNonZero = iTgtRuleCoeffs(iTgtRuleCoeffsNonZeroInd);
    
%% ROUND FOR SIGNIFICANT FIGURES
iTgtRuleCoeffsNonZeroAndSignif = round(...
    10^nSignificantFigs*iTgtRuleCoeffsNonZero)/10^nSignificantFigs;

%% CREATE COLUMN CELL STRING ARRAY WITH FIXED POINT FORMAT
% Allows for specified number of digits after decimal place.
iTgtRuleSignifCoeffStrs = strtrim(cellstr(num2str(...
    iTgtRuleCoeffsNonZeroAndSignif,...
    ['%.',num2str(nDigitsAfterDecimal),'f'])));

%% REMOVE TRAILING ZEROS
splitCell = regexp(iTgtRuleSignifCoeffStrs,'\.0+$|((?<!\.)0+$)','split');
nCoeffs = size(splitCell,1);
for iCoeff = 1:nCoeffs
    iTgtRuleSignifCoeffStrs{iCoeff} = splitCell{iCoeff}{1};
end
    
end