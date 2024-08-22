function [S,b] = create_instrument_bound_constraint_matrices(...
    instrumentMnems,ConstraintInfo)
% Helper function for OD toolkit to create matrices S and b in S*r >= b
%
% INPUTS:   
%   -> instrumentMnems: single string or nr*1 column cell string array
%   -> ConstraintInfo: sructure describing S*r >= b:
%       - instrumentMnems: single string or column cell string array
%       - instrumentCoeffs: numeric column vector of coeffs on instruments
%       - constants: numeric column vector of constraint constants
%
% OUTPUTS:  
%   -> S: nmu*nr coefficient on instruments
%   -> b: nmu*1 vector of constants
%
% This version: 05/03/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 2
    error('This function requires 2 inputs')
elseif ~ischar(instrumentMnems) && ...
        ~is_string_or_column_cell_string_array(instrumentMnems)
    error(['1st input must be string or column cell string array of ',...
        'instrument mnemonics']);
elseif ~isstruct(ConstraintInfo)
    error(['2nd input must be structure of information about the ',...
        'instrument bound constraints']);    
end

%% HANDLE OPTIONAL CHAR FORMAT OF INSTRUMENT MNEMS
if ischar(instrumentMnems)
    instrumentMnems = {instrumentMnems};
end

%% VALIDATE CONSTRAINT INFO STRUCTURE 
% Note that it allows for an optional instrumentInds fields.  This
% facilitates recycling of the ODmodel output to this function as an ODinfo
% input with some of the fields changed.
check_field_names_in_structure(...
    ConstraintInfo,'ODinfo.Constraints',...
    {'instrumentMnems';'instrumentCoeffs';'constants'},...
    {'instrumentInds';'shadowShockMnems'});

%% CHECK THE CONTENT OF THE CONSTRAINT INFO STRUCTURE
% Checks include that the dimensions of the constraint inputs are 
% consistent and that the constraint instrument mnemonics appear in the
% set of instrument mnemonics. In theory, constraint constants could be 
% -inf/inf.  In practice, this must be ruled out because 0*inf = NaN.  
constraintMnems = ConstraintInfo.instrumentMnems;
constraintCoeffs = ConstraintInfo.instrumentCoeffs;
constraintConstants = ConstraintInfo.constants;
if ~is_string_or_column_cell_string_array(constraintMnems) || ...
        ~is_finite_real_numeric_column_vector(constraintCoeffs) || ...
        ~is_finite_real_numeric_column_vector(constraintConstants)
    error(['Mnemonics of instruments to which constraints apply must ',...
        'be specified as single strings or column cell string arrays ',...
        'and constraint coefficients and constants must be specified ',...
        'as finite real and real numeric column vectors respectively'])
elseif size(constraintMnems,1)~=size(constraintCoeffs,1) || ...
        size(constraintMnems,1)~=size(constraintConstants,1)
    error(['Constraint mnemonics, coefficients and constants must be ',...
        'consistent in number (with number of rows measuring the ',...
        'number of constraints'])
elseif ~all(ismember(constraintMnems,instrumentMnems))   
    error(['At least one of the constraint instrument variable ',...
        'mnemonics does not appear to feature in the list of ',...
        'instrument mnemonics'])
end

%% LOOKUP INDICES OF CONSTRAINT INSTRUMENT MNEMONICS IN SET OF INSTRUMENTS
constraintInds = lookup_index_numbers_in_string_array(...
    instrumentMnems,constraintMnems,false);

%% CREATE COEFFICIENT MATRIX, S
nr = size(instrumentMnems,1);
nmu = size(constraintMnems,1);
S = zeros(nmu,nr);
for imu = 1:nmu
    S(imu,constraintInds(imu)) = constraintCoeffs(imu);
end

%% CREATE CONSTANT VECTOR, b
b = constraintConstants;

end