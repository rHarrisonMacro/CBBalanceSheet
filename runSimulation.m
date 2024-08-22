% This script simulates a ZLB episode in different variants of the model
% under optimal monetary policy.

%% FILE SYSTEM AND PATH MANAGEMENT
tidyUpAndSetPath;

%% START CLOCK FOR TIMING INFORMATION
clockStartTime = clock;

%% OPTIONS
variantName = 'baseline';
% Inputs
globalSolutionInputsFileName = [inputsFolder variantName '.mat'];
outputsFileName = [outputsFolder variantName '.mat'];
% Results
saveResults = true;

% Simulation options
simLength = 20;

%% LOAD RESULTS FOR THE BASELINE VARIANT
load(globalSolutionInputsFileName,'Sr','Stilde','S','pars','b');

%% SET UP EXPERIMENT SPECIFICATION
resultsFileName = [outputsFolder 'normalisationSim.mat'];
rhorstar = pars.rhorstar;
q0 = 0.75;
q0alt = q0;
load(outputsFileName,'polFunctions');
fAlt = polFunctions.notStateContingent.f;
[fSliceVarsig0,sSlice] = slice_policy_function(fAlt,Stilde,1,0);
% Slice the slice conditional on q=qHigh
fSlice = slice_policy_function(fSliceVarsig0,sSlice,2,q0);
r = fSlice(:,5);
zlbInd = find(r>b,1);
rStar0 = Sr(zlbInd-1);

%% LOAD RESULTS FOR THE BASELINE VARIANT
load(outputsFileName,'polFunctions');
f = polFunctions.stateContingent.f;
fBar = polFunctions.stateContingent.fbar;
fAlt = polFunctions.notStateContingent.f;
fAltBar = polFunctions.notStateContingent.fbar;

%% BUILD INPUTS FOR SIMULATION
%rstarSequence = Sr(rstarIndices)';
rstarSequence = rStar0*(rhorstar.^(0:simLength-1));
Ssequence = [0*rstarSequence; rstarSequence];
Snodes = get_uni_dimensional_nodes_from_tensor_product_grid(Stilde);

%% RUN SIMULATIONS 
T = length(rstarSequence);
nxBase = size(f,2);
xSimBase = nan(nxBase,T);
xSimAlt = nan(nxBase,T);
ExSimBase = nan(nxBase,T);
ExSimAlt = nan(nxBase,T);
qLag = q0;
qLagAlt = q0alt;
for t = 1:T
    % Record exogenous state vector for this period
    sSimt = Ssequence(:,t)';
    % Simulate outcomes for baseline model
    stildeSimt = [sSimt qLag];
    txBase = evaluate_endogenous_vars_in_ES_model_using_linear_interpolation(...
        Snodes,f,stildeSimt);
    EtxBase = evaluate_endogenous_vars_in_ES_model_using_linear_interpolation(...
        Snodes,fBar,stildeSimt);
    xSimBase(:,t) = txBase';
    ExSimBase(:,t) = EtxBase';
    qLag = txBase(3);
    % Simulate outcomes alternative case
    stildeSimtAlt = [sSimt qLagAlt];
    txAlt = evaluate_endogenous_vars_in_ES_model_using_linear_interpolation(...
        Snodes,fAlt,stildeSimtAlt);
    EtxAlt = evaluate_endogenous_vars_in_ES_model_using_linear_interpolation(...
        Snodes,fAltBar,stildeSimtAlt);
    xSimAlt(:,t) = txAlt';
    ExSimAlt(:,t) = EtxAlt';
    qLagAlt = txAlt(3);
end

%% SAVE RESULTS
if saveResults
    save(resultsFileName,'xSimBase','xSimAlt','ExSimBase','ExSimAlt',...
        'rstarSequence','pars','q0','q0alt');
end

%% REPORT TIME TAKEN
clockFinishTime = clock;
compute_and_display_elapsed_time(clockStartTime,clockFinishTime);