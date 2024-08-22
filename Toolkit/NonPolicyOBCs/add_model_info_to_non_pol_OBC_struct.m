function nonPolOBCs = add_model_info_to_non_pol_OBC_struct(nonPolOBCs,...
    Model,OptPolInfo)
% This function adds model-specific information about the effects of
% non-instrument constraints to the structure describing those constraints.
% This information is used by lower-level functions to construct the
% time-varying matrices used to implement the algorithms in Section 6 of
% "Optimal policy with occasionally binding constraints: piecewise linear
% solution methods" by Richard Harrison and Matt Waldron.
% 
% INPUTS
%   - nonPolOBCs: structure containing high level information about the
%     non-instrument constraints. These are in fields "active" and
%     "inactive", which both contain a field "VarMnems", containing the
%     mnemnomics of the variables that are constrained when constraints are
%     active or inactive. The field "active" also contains a "Constants"
%     field, which contains the values that the constrained variables must
%     satisfy when the constraints are active.
%
% OUTPUTS
%   - nonPolOBCs: structure containing high level information about the
%     non-instrument constraints, augmented to included lower level
%     information:
%       - EqInds: indices of the equations that define when each constraint 
%         is "inactive".
%       - G: N*nxtilde matrix indicating which variable is bounded when
%         each constraints is "active".
%       - k: N*1 vector of bound values for constrained variables
%       - active.varInds: indices of variables subject to bounds when
%       constraints are active.
%       - inactive.varInds: indices of variables subject to bounds when
%       constraints are inactive.
%       - xTildeInds: indices of non-policy variables (xtilde) within full
%       set of model variables (x).
%
% This version: 10/12/2020
% Author: Richard Harrison


%% PARTITION THE MODEL STRUCTURE ACCORDING TO OPTIMAL POLICY INFORMATION
[HB,HC,HF,PSI,xMnems,zMnems,xEqNames] = unpack_model(Model,...
    {'HB';'HC';'HF';'PSI';'xMnems';'zMnems';'xEqNames'});

[HtildeBxtilde,~,~,HtildeFxtilde,~,PSItildeztilde,instrumentLogicals,...
    ~,policyEqLogicals] = ...
    partition_structural_equations_for_optimal_discretion_solution(...
    HB,HC,HF,PSI,xMnems,zMnems,xEqNames,OptPolInfo.instrumentMnems,...
    OptPolInfo.policyEqNames,OptPolInfo.policyShockMnems,...
    OptPolInfo.objVarMnems);

%% ADD INDEXING INFORMATION TO NON-POLICY OBC INFO
N = size(nonPolOBCs.inactive.VarMnems,1);
xtildeEqNames = xEqNames(~policyEqLogicals);
nxtilde = size(xtildeEqNames,1);
xtildeMnems = xMnems(~instrumentLogicals);
xtildeInds = 1:size(xMnems,1);
xtildeInds(instrumentLogicals) = [];
PSIdelta = zeros(nxtilde,N);
EqInds = nan(N,1);
G = zeros(N,nxtilde);
activeVarInds = nan(N,1);
inactiveVarInds = nan(N,1);
for iOBC = 1:N
    iEqName = nonPolOBCs.inactive.EqNames{iOBC};
    iEqInd = lookup_model_index_numbers(xtildeEqNames,iEqName);
    PSIdelta(iEqInd,iOBC) = nonPolOBCs.active.Constants(iOBC);
    EqInds(iOBC) = iEqInd;
    iActiveVarMnem = nonPolOBCs.active.VarMnems{iOBC};
    iActiveVarInd = lookup_model_index_numbers(xtildeMnems,iActiveVarMnem);
    activeVarInds(iOBC) = iActiveVarInd;
    G(iOBC,iActiveVarInd) = 1;
    iInactiveVarMnem = nonPolOBCs.inactive.VarMnems{iOBC};
    inactiveVarInds(iOBC) = ...
        lookup_model_index_numbers(xtildeMnems,iInactiveVarMnem);
end

%% LOAD RESULTS BACK INTO STRUCTURE
nonPolOBCs.EqInds = EqInds;
nonPolOBCs.PSIdelta = PSIdelta;
nonPolOBCs.G = G;
nonPolOBCs.k = nonPolOBCs.active.Constants;
nonPolOBCs.active.VarInds = activeVarInds;
nonPolOBCs.inactive.VarInds = inactiveVarInds;
nonPolOBCs.xtildeInds = xtildeInds;

%% PERFORM VALIDATION CHECKS ON THE EQUATIONS THAT SHOULD BE DELETED
if any(any(HtildeBxtilde(EqInds,:)~=0)) ||  ...
        any(any(HtildeFxtilde(EqInds,:)~=0)) || ...
        any(any(PSItildeztilde(EqInds,:)~=0))
    error(['Equations to be removed from the system when non-policy ',...
        'OBCs are binding cannot contain leads, lags or loadings on ',...
        'shocks.']);
end

