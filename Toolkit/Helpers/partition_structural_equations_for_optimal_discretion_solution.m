function [HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,instrumentLogicals,instrumentInds,policyEqLogicals,...
    policyEqInds,policyShockLogicals,policyShockInds] = ...
    partition_structural_equations_for_optimal_discretion_solution(...
    HB,HC,HF,PSI,xMnems,zMnems,xEqNames,instrumentMnems,policyEqNames,...
    policyShockMnems,objVarMnems)
% Helper function for OD toolkit to partition model equations.
%
% INPUTS:   
%   -> HB: backward loadings in structural equations
%   -> HC: contemporaneous loadings in structural equations
%   -> HF: forward loadings in structural equations
%   -> PSI: shock loadings in structural equations
%   -> xMnems: model variable mnemonics
%   -> zMnems: shock mnemonics
%   -> xEqNames: equation names
%   -> instrumentMnems: instrument mnemonics
%   -> policyEqNames: instrument equation names
%   -> policyShockMnems: instrument shock mnemonics
%   -> objVarMnems: single string or column cell string array
%
% OUTPUTS:  
%   -> HtildeBxtilde: HB with policy rule row and instrument(s) col removed
%   -> HtildeCxtilde: HC with policy rule row and instrument(s) col removed
%   -> HtildeCr: HC with policy rule row and non instrument(s) col removed
%   -> HtildeFxtilde: HF with policy rule row and instrument(s) col removed
%   -> HtildeFr: HF with policy rule row and non instrument(s) col removed
%   -> PSItildeztilde: PSI with policy rule row and shock col removed
%       - policyEqLogicals: vector of logicals that locate policy eqs
%       - policyEqInds: vector of indices that locate policy eqs
%       - policyShockLogicals: vector of logicals locating policy shocks
%       - policyShockInds: vector of indices locating policy shocks
%       - instrumentLogicals: vector of logicals locating the instruments
%       - instrumentInds: vector of indices that locate the instruments
%
% This version: 06/03/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
% The checks include that the input optimisation information is consistent
% with each other and with the model being used.
if nargin < 11
    error('This function requires 11 inputs')
elseif ~is_string_or_column_cell_string_array(instrumentMnems) || ...
       ~is_string_or_column_cell_string_array(policyEqNames) || ... 
       ~is_string_or_column_cell_string_array(policyShockMnems) || ... 
       ~is_string_or_column_cell_string_array(objVarMnems) 
    error(['Objective mnemonics, instrument mnemonics, policy ',...
        'shock mnemonics and policy equation names must be ',...
        'specified as single strings or column cell string arrays'])
elseif size(instrumentMnems,1)~=size(policyEqNames,1)
    error(['You must specify the same number of instrument mnemonics ',...
        'and policy equations'])
elseif ~all(ismember(objVarMnems,xMnems)) || ...
        ~all(ismember(instrumentMnems,xMnems))    
    error(['At least one objective variable or instrument variable ',...
        'mnemonic does not appear to feature in the model being used'])
elseif ~all(ismember(policyEqNames,xEqNames)) 
    error(['At least one policy equation name does not appear to ',...
        'feature in the model being used'])  
elseif ~isempty(policyShockMnems) && ...
        ~all(ismember(policyShockMnems,zMnems))
    error(['At least one policy shock mnemonic does not appear to ',...
        'feature in the model being used'])      
elseif any(ismember(objVarMnems,instrumentMnems))
    error('The list of objective and instrument mnemonics overlaps');
end

%% COMPUTE INSTRUMENT, POLICY EQUATION AND POLICY SHOCK LOGICALS
instrumentLogicals = ismember(xMnems,instrumentMnems);
policyEqLogicals = ismember(xEqNames,policyEqNames);
policyShockLogicals = ismember(zMnems,policyShockMnems);

%% COMPUTE INSTRUMENT, POLICY EQUATION AND POLICY SHOCK INDICES
instrumentInds = lookup_index_numbers_in_string_array(...
    xMnems,instrumentMnems);
policyEqInds = lookup_index_numbers_in_string_array(...
    xEqNames,policyEqNames);
if any(policyShockLogicals)
    policyShockInds = lookup_index_numbers_in_string_array(...
        zMnems,policyShockMnems);
else
    policyShockInds = [];
end

%% PARTITION THE STRUCTURAL EQUATIONS
% Note the use of instrument indices to take care of the ordering of the
% instrument mnemonics.
HtildeBxtilde = HB(~policyEqLogicals,~instrumentLogicals);
HtildeBr = HB(~policyEqLogicals,instrumentInds);
HtildeCxtilde = HC(~policyEqLogicals,~instrumentLogicals);
HtildeCr = HC(~policyEqLogicals,instrumentInds);
HtildeFxtilde = HF(~policyEqLogicals,~instrumentLogicals);
HtildeFr = HF(~policyEqLogicals,instrumentInds);
PSItildeztilde = PSI(~policyEqLogicals,~policyShockLogicals);

%% CHECK THAT THERE ARE NO LAGGED INSTRUMENTS IN THE NON-INSTRUMENT EQS
% This is an assumption in the MAPS algorithm and the result is not valid
% if this assumption is violated.
if any(any(abs(HtildeBr)>1e-12))
    error(['One or more lags of the instruments appear in the model ',...
        'equations.  The MAPS algorithms assume and requires that ',...
        'this is not the case.  There are two ways of fixing this ',...
        'problem: i) create a new variable with an identity equation ',...
        'that you can then use as the instrument; ii) swap one of the ',...
        'objective variables with the instrument'])
end

end