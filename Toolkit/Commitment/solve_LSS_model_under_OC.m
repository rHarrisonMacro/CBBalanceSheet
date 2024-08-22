function [By,PHIy,Fy,HBy,HCy,HFy,PSIy,OCinfo] = ...
    solve_LSS_model_under_OC(Model,OCinfo) 
% This function solves a RE model under commitment using the approach of
% Dennis (2007) "Optimal policy in rational expectations models: new
% solution algorithms", Macroeconomic Dynamics, 11(1).
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
%       - Constraints (optional): structure containing information about
%       constraints on the policy instruments
%
% OUTPUTS:  
%   -> BOC: transition for endogenous variables
%   -> PHIOC: shock impact
%   -> FOC: forward loadings
%   -> HBOC: coeeficients on lags in structural form
%   -> HCOC: coefficients on contemporaneous variables in structural form
%   -> HFOC: coefficients on leads in structural form
%   -> PSIOC: coefficients on shocks in structural form
%   -> OCinfo: input structure updated with:
%       - policyEqLogicals: vector of logicals that locate policy eqs
%       - policyEqInds: vector of indices that locate policy eqs
%       - policyShockLogicals: vector of logicals locating policy shocks
%       - policyShockInds: vector of indices locating policy shocks
%       - instrumentLogicals: vector of logicals locating the instruments
%       - instrumentInds: vector of indices that locate the instruments
%       - Q: instruments weights matrix
%       - W: target weights matrix
%
%
% This version: 25/07/2019
% Author(s): Matt Waldron & Richard Harrison

%% CHECK INPUTS
if nargin < 2
    error('This function requires 2 or 3 inputs')
end

%% UNPACK AND VALIDATE OPTIONS INPUT
% TBC

%% UNPACK MODEL
[HB,HC,HF,PSI,xMnems,zMnems,xEqNames] = unpack_model(...
    Model,{'HB','HC','HF','PSI','xMnems','zMnems','xEqNames'});

%% UNPACK THE INFO STRUCTURE
policyEqNames = OCinfo.policyEqNames;
policyShockMnems = OCinfo.policyShockMnems;
instrumentMnems = OCinfo.instrumentMnems;
instrumentWeights = OCinfo.instrumentWeights;
objVarMnems = OCinfo.objVarMnems;
objVarWeights = OCinfo.objVarWeights;
beta = OCinfo.beta;

%% PARTITION THE STRUCTURAL EQUATIONS
[HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,instrumentLogicals,instrumentInds,policyEqLogicals,...
    policyEqInds,policyShockLogicals,policyShockInds] = ...
    partition_structural_equations_for_optimal_discretion_solution(...
    HB,HC,HF,PSI,xMnems,zMnems,xEqNames,instrumentMnems,policyEqNames,...
    policyShockMnems,objVarMnems);

%% CREATE ENDOGENOUS VARIABLES AND INSTRUMENTS LOSS FUNCTION WEIGHTS
[W,Q] = create_loss_function_weight_matrices(...
    xMnems,objVarMnems,objVarWeights,instrumentMnems,instrumentWeights);

%% UPDATE THE OCINFO STRUCTURE WITH THE LOGICALS AND LOSS FUNCTION WEIGHTS
OCinfo.instrumentLogicals = instrumentLogicals;
OCinfo.instrumentInds = instrumentInds;
OCinfo.policyEqLogicals = policyEqLogicals;
OCinfo.policyEqInds = policyEqInds;
OCinfo.policyShockLogicals = policyShockLogicals;
OCinfo.policyShockInds = policyShockInds;
OCinfo.Q = Q;
OCinfo.W = W;

%% ADJUST SHOCK LOADING MATRICES TO INCLUDE ZERO LOADINGS ON POLICY SHOCKS
nz = size(PSI,2);
nxTilde = size(HtildeBxtilde,1);
PSItildeztildeWithPolShocks = zeros(nxTilde,nz);
PSItildeztildeWithPolShocks(:,~policyShockLogicals) = PSItildeztilde;

%% COMPUTE SOLUTION USING DENNIS ALGORITHM REFACTORED TO MAPS NOTATION
[By,PHIy,Fy,HBy,HCy,HFy,PSIy,OCinfo] = ...
    compute_OC_solution_using_Dennis_algorithm(...
    HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztildeWithPolShocks,beta,OCinfo);


end