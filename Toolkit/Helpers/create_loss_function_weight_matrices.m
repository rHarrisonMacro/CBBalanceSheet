function [W,Q] = create_loss_function_weight_matrices(...
    xMnems,objVarMnems,objVarWeights,instrumentMnems,instrumentWeights)
% Helper function for OD toolkit to create weights for loss function.
%
% INPUTS:   
%   -> xMnems: model variable mnemonics
%   -> objVarMnems: single string or column cell string array
%   -> objVarWeights: column vector of loss weights for objective variables
%   -> instrumentMnems: single string or column cell string array
%   -> instrumentWeights: column vector of loss weights for instruments
%
% OUTPUTS:  
%   -> W: weights on endogenous variables in loss
%   -> Q: weights on instruments in loss
%
% This version: 06/03/2017
% Author(s): Matt Waldron

%% CHECK INPUTS
% Checks include that the dimensions of the instrument and objective
% metadata and weights are consistent, that the metadata appears in the
% model being used, and that there is no overlap between the instrument and
% objective mnemonic groups.
if nargin < 5
    error('This function requires 5 inputs')
elseif ~is_string_or_column_cell_string_array(objVarMnems) || ...
        ~is_finite_non_negative_numeric_column_vector(objVarWeights) || ...
        ~is_string_or_column_cell_string_array(instrumentMnems) || ...
        ~is_finite_non_negative_numeric_column_vector(instrumentWeights)
    error(['Objective and instrument mnemonics must be specified as ',...
        'single strings or column cell string arrays and objective ',...
        'and instrument weights must be specified as non-negative ',...
        'finite numeric column vectors'])
elseif size(objVarMnems,1)~=size(objVarWeights,1) || ...
        size(instrumentMnems,1)~=size(instrumentWeights,1)
    error(['Objective and instrument mnemonics and weights must be ',...
        'consistent in number'])
elseif ~all(ismember(objVarMnems,xMnems)) || ...
        ~all(ismember(instrumentMnems,xMnems))    
    error(['At least one objective variable or instrument variable ',...
        'mnemonic does not appear to feature in the model being used'])
elseif any(ismember(objVarMnems,instrumentMnems))
    error('The list of objective and instrument mnemonics overlaps');
end

%% CREATE Q
Q = diag(instrumentWeights);

%% COMPUTE DIMENSIONALITY OF NON-INSTRUMENT ENDOGENOUS VARIABLES
nx = size(xMnems,1);
nr = size(instrumentMnems,1);
nxtilde = nx-nr;

%% COMPUTE INSTRUMENT LOGICALS
instrumentLogicals = ismember(xMnems,instrumentMnems);

%% CREATE LIST OF NON-INSTRUMENT ENDOGENOUS VARIABLE MNEMONICS
xtildeMnems = xMnems(~instrumentLogicals);

%% LOOKUP INDEX NUMBERS OF OBJECTIVE VARS IN NON-INSTRUMENT ENDOG VARS
xtildeInds = lookup_model_index_numbers(xtildeMnems,objVarMnems);

%% CREATE W
W = zeros(nxtilde,nxtilde);
W(xtildeInds,xtildeInds) = diag(objVarWeights);

end