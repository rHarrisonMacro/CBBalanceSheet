function [f0,Ef0] = ...
    compute_piecewise_linear_solution_block(ModelXiL,...
        ModelXiH,OptPolInfo,Stilde)

% This code is computes the piecewise linear solution as an initial guess 
% for the non-linear policy function.
%
%% SET INITIAL OPTIONS
[nSims,nstildePlus1] = size(Stilde);
% nz = nstildePlus1 - 1;
simHorizon = 25;

ODPPopts = struct;
ODPPopts.Projection.BPY.useBestOnNonConvergence = true;
ODPPopts.Projection.BPY.issueCyclingWarning = false;

ODPPopts.Solution.Algorithm.maxIter = 20000;
ODPPopts.Solution.Algorithm.tol = 1e-06;

% BPYopts = struct;
% % TBC set up options to get best performance
% BPYopts.useBestOnNonConvergence = true;
% BPYopts.issueCyclingWarning = false;

endogVarMnems = {'pie';'x';'q';'mux';'ir';'muxL';'muxH';'varsigChk'};
nEndogVars = size(endogVarMnems,1);

exogStateVarMnems = {'varsigma';'rstar'};
nExogStates = size(exogStateVarMnems,1);
exogStateShocks = cell(nExogStates,1);
exogStateSigmas = cell(nExogStates,1);
for iState = 1:nExogStates
    exogStateShocks{iState} = ['eta' exogStateVarMnems{iState}];
    exogStateSigmas{iState} = ['sigma' exogStateVarMnems{iState}];
end

%% UNPACK MODEL INFO
[theta,thetaMnems,xMnems,zMnems] = ...
    unpack_model(ModelXiL,{'theta';'thetaMnems';'xMnems';'zMnems'});
endogVarInds = lookup_model_index_numbers(xMnems,endogVarMnems);
exogStateShkInds = lookup_model_index_numbers(zMnems,exogStateShocks);

exogStateSigmaInds = ...
    lookup_model_index_numbers(thetaMnems,exogStateSigmas);
exogStateSigmaVals = theta(exogStateSigmaInds);

nz = size(zMnems,1);
nx = size(xMnems,1);

q0vec = Stilde(:,end);
qInd = lookup_model_index_numbers(xMnems,'q');

% %% GENERATE BASIC INPUTS FOR PROJECTION CODES BELOW
% % Solve for time-invariant OD solution (using default options)
% [BOD,PHIOD,HBOD,HCOD,HFOD,PSIOD,~,OptPolInfo,Vxtildextilde] = ...
%     solve_LSS_model_under_optimal_discretion(Model,OptPolInfo);
% % Extract info from OptPolInfo structure
% W = OptPolInfo.W;
% Q = OptPolInfo.Q;
% beta = OptPolInfo.beta;
% instrumentLogicals = OptPolInfo.instrumentLogicals;
% policyEqLogicals = OptPolInfo.policyEqLogicals;
% instrumentInds = OptPolInfo.instrumentInds;
% policyShockLogicals = OptPolInfo.policyShockLogicals;
% % Extract partitioned matrices for BPY
% Bxtildextilde = BOD(~instrumentLogicals,~instrumentLogicals);
% Brxtilde = BOD(instrumentInds,~instrumentLogicals);
% HtildeBxtilde = HBOD(~policyEqLogicals,~instrumentLogicals);
% HtildeCxtilde = HCOD(~policyEqLogicals,~instrumentLogicals);
% HtildeCr = HCOD(~policyEqLogicals,instrumentInds);
% HtildeFxtilde = HFOD(~policyEqLogicals,~instrumentLogicals);
% HtildeFr = HFOD(~policyEqLogicals,instrumentInds);
% PSItildeztilde = PSIOD(~policyEqLogicals,~policyShockLogicals);
% % Create a dummy F matrix (since there are no anticipated shocks)
% FODdummy = 0*BOD;
% % Extract information about instrument bounds
% instrumentMnems = OptPolInfo.instrumentMnems;
% instrumentInds = ...
%     lookup_model_index_numbers(xMnems,instrumentMnems);
% [S,b] = create_instrument_bound_constraint_matrices(...
%     instrumentMnems,OptPolInfo.Constraints);
% nmu = size(b,1);
% bmat = repmat(b,1,simHorizon);
% 
% % Compute dimensions of shocks and endogenous variables
% [nx,nz] = size(PHIOD);
% 
% Create a baseline set of shocks
tShocks = struct;
tShocks.anticipated = zeros(nz,simHorizon);

% Prepare for loop
f0 = nan(nSims,nEndogVars);
Ef0 = nan(nSims,nEndogVars);
x0 = zeros(nx,1);

for t = 1:nSims
    jS = Stilde(t,1:2)';
    tShocks.anticipated(exogStateShkInds,1) = jS./exogStateSigmaVals;
    x0(qInd) = q0vec(t);
    if jS(1)==0
        ix = compute_ODPP(ModelXiL,OptPolInfo,x0,tShocks,ODPPopts);
    else
        ix = compute_ODPP(ModelXiH,OptPolInfo,x0,tShocks,ODPPopts);
    end
    f0(t,:) = ix(endogVarInds,1)';
    Ef0(t,:) = ix(endogVarInds,2)';
end


end

