function [xfeConstrained,areConstraintsBinding,mufeh] = ...
    project_LSS_model_variables_subject_to_constraints(Model,x0,Shocks,...
    ConstrainedSimulInfo)
% Constructs an LSS model projection with occasionally binding constraints,
% implemented using the method of Holden and Paetz (2012) "Efficient
% simulation of DSGE models with inequality constraints".
% 
% INPUTS
%   -> Model: MAPS LSS model structure
%   -> x0: nx*1 vector of initial conditions
%   -> Shocks: structure containing the following fields:
%       - anticipated: nz*H matrix of anticipated shock values
%   -> ConstrainedSimulInfo: structure containing the following fields
%       - instrumentVarMnems: mnemonics of variables subject to constraints
%       - shadowShockMnems: mnemonics of shocks used to implement
%       constraints.
%       - boundVals: nmu*1 vector of bound values for the constraints
%       - expectationsHorizon (optional): integer scalar for number of 
%       additional projection periods to use to apply bounds.
%       - HPopts (optional). Structure containing options for Holden-Paetz
%       solution:
%           - mu0 (optional): initial condition for the search (must be 
%           specified as an nmu*H matrix (rows indexing shadow shock, 
%           cols periods)
%           - fvalTol (optional): tolerance for accepting the Linear 
%           Complementarity Solution
%           - Quadprog (optional): structure of options for the algorithm
%           (see invert_to_find_shadow_shocks_using_Holden_Paetz.m)

%
% OUTPUTS
%   -> xfeConstrained: nx*H projection of endogenous variables, subject to
%      the occasionally binding constraints.
%   -> areConstraintsBinding: scalar logical indicating whether or not any
%      constraints bind in the simulation.
%
%
% This version: 20/05/2019
% Authors: Richard Harrison and Matt Waldron

%% UNPACK CONSTRAINED SIMULATION INFO (IF APPLICABLE)
[xmuInds,zmuInds,b,Htail,HPopts] = ...
    check_and_unpack_constraint_info_structure(...
    Model,ConstrainedSimulInfo);

%% CHECK SHOCK STRUCTURE
if isfield(Shocks,'unanticipated')
    error(['This function cannot be used with unanticipated shocks, ',...
        'except in the first period of projection (which can be ',...
        'incorporated in the first column of the anticipated shocks ',...
        'matrix.']);
end
if ~isfield(Shocks,'anticipated')
    error(['This function requires an anticipated shocks matrix within ',...
        'the ''Shocks'' structure.']);
end

%% GATHER INFORMATION FROM INPUTS
H = size(Shocks.anticipated,2);
Hfe = H + Htail;

%% UNPACK LSS MODEL INFO
[B,PHI,F] = unpack_model(Model,{'B','PHI','F'});

%% ANALYSE DIMENSIONS OF PROBLEM
nz = size(PHI,2);

%% EXTRACT SHOCKS TO APPLY
ShocksUnconstrained = struct;
ShocksUnconstrained.anticipated = [Shocks.anticipated, zeros(nz,Htail)];

%% COMPUTE AN UNCONSTRAINED PROJECTION
xfeUnconstrained = ...
    project_LSS_model_variables(B,PHI,F,x0,ShocksUnconstrained);
xmufeUnconstrained = xfeUnconstrained(xmuInds,:);

%% CHECK IF CONSTRAINTS BIND
% This repeats the logic of a check in the lower-level function called
% below. However, implementing this check at this point allows the code to
% test whether constraints are binding before executing the potentially
% expensive step of constructing the multiplier matrices.
T = size(ShocksUnconstrained.anticipated,2);
xfminusBounds = xmufeUnconstrained'-ones(T,1)*b';
v = xfminusBounds(:);
if all(v>=0)
    areConstraintsBinding = false;
    xfeConstrained = xfeUnconstrained;
    nmu = size(zmuInds,1);
    mufeh = zeros(nmu,Hfe);
    return
end

%% CONSTRUCT MULTIPLIER MATRIX FOR H-P INVERSION
M = construct_shadow_shocks_impact_matrix_for_Holden_Paetz(...
    B,PHI,F,xmuInds,zmuInds,Hfe);

%% COMPUTE HOLDEN-PAETZ INVERSION
[mufeh,areConstraintsBinding] = ...
    invert_to_find_shadow_shocks_using_Holden_Paetz(...
            xmufeUnconstrained,b,M,HPopts);

%% COMPUTE CONSTRAINED PROJECTION
if areConstraintsBinding
    ShocksConstrained = ShocksUnconstrained;
    ShocksConstrained.anticipated(zmuInds,:) = mufeh;
    xfeConstrained = project_LSS_model_variables(...
        B,PHI,F,x0,ShocksConstrained);
else
    xfeConstrained = xfeUnconstrained;
end
       
end

%% HELPER FUNCTION TO HANDLE CONSTRAINT INFO STRUCTURE
function [xmuInds,zmuInds,b,Htail,HPopts] = ...
    check_and_unpack_constraint_info_structure(Model,ConstrainedSimulInfo)
% Computes a stochastic simulation with LSS model.
%
% INPUTS:
%   -> Model: LSS model
%   -> ConstrainedSimulInfo:
%       - instrumentVarMnems: mnemonics of instruments subject to bound 
%         constraints
%       - shadowShockMnems: mnemonics of shadow shocks to impose those
%         constraints
%       - boundVals: nBounds*1 vector of bounds for the instruments
%       - expectationsHorizon: number of periods over which to compute 
%         expectations
%       - HPopts (optional): options for the HP algorithm
%           - fvalTol: tolerance for convergence
%
% OUTPUTS
%   -> xmuInds: indices of instruments subject to bound constraints
%   -> zmuInds: indices of shadow shocks to use to impose the,
%   -> b: nBounds*1 vector of bounds for the instruments
%   -> Htail: number of periods over which to compute expectations
%   -> HPopts: options for the HP algorithm
%       - fvalTol: tolerance for convergence

%% CHECK FIELD NAMES
compulsoryFields = {...
    'instrumentVarMnems';
    'shadowShockMnems';
    'boundVals'};
optionalFields = {...
    'HPopts'
    'expectationsHorizon'};
structureDescription = 'Constrained simulation info';
check_field_names_in_structure(ConstrainedSimulInfo,...
    structureDescription,compulsoryFields,optionalFields);

%% UNPACK STRUCTURE
instrumentVarMnems = ConstrainedSimulInfo.instrumentVarMnems;
shadowShockMnems = ConstrainedSimulInfo.shadowShockMnems;
b = ConstrainedSimulInfo.boundVals;
if isfield(ConstrainedSimulInfo,'HPopts')
    HPopts = ConstrainedSimulInfo.HPopts;
else
    HPopts = struct;
end
if isfield(ConstrainedSimulInfo,'expectationsHorizon')
    Htail = ConstrainedSimulInfo.expectationsHorizon;
else
    Htail = 0;
end

%% COMPUTE INDICES OF INSTRUMENTS
xMnems = unpack_model(Model,'xMnems');
xmuInds = lookup_model_index_numbers(xMnems,instrumentVarMnems);

%% COMPUTE INDICES OF SHADOW SHOCKS
zMnems = unpack_model(Model,'zMnems');
zmuInds = lookup_model_index_numbers(zMnems,shadowShockMnems);

%% CHECK BOUNDS
if ~is_finite_real_numeric_column_vector(b)
    error(['Constraint bound values must be represented as a finite ',...
        'real column vector']);
end

%% CHECK CONSISTENCY OF BOUND INFO
if size(xmuInds,1)~=size(zmuInds,1) || size(xmuInds,1)~=size(b,1)
    error(['Number of instruments subject to bound constraints, ',...
        'shadow shocks used to implement bound constraints and bound ',...
        'constraint values must be equal']);
end

%% CHECK EXPECTATIONS HORIZON
if ~is_non_negative_real_integer(Htail)
    error('Expectations horizon must be a finite non-negative integer.');
end
    
end

