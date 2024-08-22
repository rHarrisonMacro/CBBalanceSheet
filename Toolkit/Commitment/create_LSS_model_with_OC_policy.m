function [ModelwOCpolicy,OCinfo] = ...
    create_LSS_model_with_OC_policy(Model,OCinfo)
% This function creates a temporary MAPS model (which is *not* a valid MAPS
% model structure) constaining numeric matrices for the structural form and
% RE solution for the policy problem augmented with Lagrange multipliers.
%
% It is a wrapper for solve_LSS_model_under_OC(Model,OCinfo)
%
%
% INPUTS:   
%   -> Model: MAPS LSS model
%   -> OCinfo: structure of info required to solve the OC problem
%       - policyEqNames: single string or column cell string array
%       - policyShockMnems: single string or column cell string array (use
%         empty string ('') if there are no policy shocks to remove)
%       - instrumentMnems: single string or column cell string array
%       - instrumentWeights: column vector of loss weights for instruments
%       - objVarMnems: single string or column cell string array
%       - objVarWeights: column vector of loss weights for objective vars
%       - beta: discount factor
%       - Constraints: (optional) structure containing information about
%       instrument constraints (bounds)
%
% OUTPUTS:  
%   -> ModelwOCpolicy: model with optimal commitment policy -- nb. this is
%      not a valid MAPS model
%
%
% This version: 25/07/2019
% Author(s): Matt Waldron & Richard Harrison


%% SOLVE MODEL
[By,PHIy,Fy,HBy,HCy,HFy,PSIy,OCinfo] = ...
    solve_LSS_model_under_OC(Model,OCinfo);

%% AUGMENT METADATA
xMnems = unpack_model(Model,'xMnems');
xtildeMnems = xMnems(~OCinfo.instrumentLogicals);
nxtilde = size(xtildeMnems,1);
lambdaMnems = cell(nxtilde,1);
for ixtilde = 1:nxtilde
    lambdaMnems{ixtilde} = ['Lagrange multiplier ',num2str(ixtilde)];
end
xMnemsOC = [xMnems;lambdaMnems];

%% ADJUST SHOCK METADATA IF NECESSARY
zMnems = unpack_model(Model,'zMnems');
if isfield(OCinfo,'Constraints')
   zMnemsOC = [zMnems;OCinfo.Constraints.shadowShockMnems];
else
   zMnemsOC = zMnems;
end

%% PACK MATRICES INTO INPUT MODEL
ModelwOCpolicy = Model;
ModelwOCpolicy = pack_model(ModelwOCpolicy,...
    {'B','PHI','F','HB','HC','HF','PSI','xMnems','zMnems'},...
    {By,PHIy,Fy,HBy,HCy,HFy,PSIy,xMnemsOC,zMnemsOC});

end