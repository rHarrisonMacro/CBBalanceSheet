% This script sets up the inputs for the global solution, starting from the
% MAPS model representation.

%% WORKSPACE AND PATH MANAGEMENT
tidyUpAndSetPath;

%% START CLOCK FOR TIMING INFORMATION
clockStartTime = clock;

%% OPTIONS
saveInputs = true;
saveFileFolder = inputsFolder;
saveFileName = 'baseline'; 

%% OPTIONS PERTAINING TO COMPUTATION OF INITIAL POLICY FUNCTION GUESS 
% (USING PIECEWISE LINEAR SOLUTION)
nBlocksForGuess = 32;

%% LOWER LEVEL OPTIONS
checkConditionalOutcomes = false;

%% INCORPORATE STATE CONTINGENCY INFO
% "Crisis" shock
% There are two states:
% -- In the "normal" state, spreads (varsigma) are zero.
% -- In the "crisis" state, spreads are high. 
% The probability of N->C transition is p. The probability of C->C
% transition is m.
p = 0.01;
m = 0.5;
varsigmaHigh = 1.05^0.25-1; % 5% annualised spread
nu = 0.0015;
xiStateContingencyFactor = 5; % Multiplicative flow QE factor (=xiH/xiL)
xiL = 0.065/xiStateContingencyFactor; 

xiH = xiL*xiStateContingencyFactor; 

%% PARAMETER ASSUMPTIONS
qLB = 0;
qUB = 0.75; 
% Weight on QE and change in loss function:
epsilon = 0.01/100; 
epsilonDelta = epsilon*20;

%% OPTIONS FOR DIAGNOSTIC PLOTS
qVal1 = 0;
qVal2 = qUB;

%% SPECIFY INFORMATION ABOUT MARKOV STATES
nr = 51; 
nq = 100; 
nvarsig = 2;
Svarsig = [0; varsigmaHigh];
Omegavarsig = [1-p, p; 1-m, m];

%% SPECIFY MODEL
MAPSmodelFileName = 'QEmodel.maps';

%% SPECIFY INFORMATION ABOUT VARIABLES IN GLOBAL SOLUTION
endogVarMnems = {'pie';'x';'q';'mux';'ir';'muxL';'muxH';'varsigChk'};
nEndogVars = size(endogVarMnems,1);
endogStateVarMnems = {'q'};
exogStateVarMnems = {'varsigma';'rstar'};
nExogStates = size(exogStateVarMnems,1);
exogStateShocks = cell(nExogStates,1);
exogStateSigmas = cell(nExogStates,1);
for iState = 1:nExogStates
    exogStateShocks{iState} = ['eta' exogStateVarMnems{iState}];
    exogStateSigmas{iState} = ['sigma' exogStateVarMnems{iState}];
end

%% SPECIFY HORIZON FOR PIECEWISE LINEAR GUESS CONSTRUCTION
H = 25;
ODPPopts = struct;
ODPPopts.Projection.BPY.useBestOnNonConvergence = true;

%% LOAD MODEL AND UPDATE PARAMETER VALUES IF REQUIRED
Model = create_model(MAPSmodelFileName);

% Re-solve two model variants depending on xiH and xiL:
ModelXiL = update_parameters_and_resolve_LSS_model(Model,...
    [0;xiL;nu],{'rhovarsigma';'xi';'nu'});

ModelXiH = update_parameters_and_resolve_LSS_model(Model,...
    [m;xiH;nu],{'rhovarsigma';'xi';'nu'});

[theta,thetaMnems,xMnems,zMnems] = ...
    unpack_model(Model,{'theta';'thetaMnems';'xMnems';'zMnems'});
endogVarInds = lookup_model_index_numbers(xMnems,endogVarMnems);
exogStateShkInds = lookup_model_index_numbers(zMnems,exogStateShocks);
endogStateVarInds = lookup_model_index_numbers(xMnems,endogStateVarMnems);
exogStateSigmaInds = ...
    lookup_model_index_numbers(thetaMnems,exogStateSigmas);
exogStateSigmaVals = theta(exogStateSigmaInds);
nx = size(xMnems,1);
nz = size(zMnems,1);

%% BUILD AUXILLARY PARAMETERS
pars = unpack_model_parameter_values_as_struct(Model);
omegax = pars.omegax;
omegapi = pars.omegapi;
omegaq = epsilon;
omegaDq = epsilonDelta;
b = log(pars.beta);

%% BUILD MARKOV PROCESS INPUTS
[Sr,Omegar] = rouwenhorst2(nr,0,pars.rhorstar,pars.sigmarstar);
Omega = kron(Omegar,Omegavarsig);
S = grid2d({Svarsig,Sr});
Sq = linspace(qLB,qUB,nq)';
Stilde = grid2d({Svarsig,Sr,Sq});
ns = nvarsig*nr;
nstilde = nq*ns;

%% BUILD ENDOGENOUS STATE LOOKUP MATRIX, LAMBDA
oneTonq = (1:1:nq);
jIndices = repmat(oneTonq,ns,1);
jIndices = jIndices(:);
Lambda = nan(nstilde,nq);
for k = 1:nstilde
    j = jIndices(k);
    for m = 1:nq
        Lambda(k,m) = k - (j-m)*ns;
    end
end

%% CONSTRUCT INFORMATION REQUIED TO COMPUTE GUESSES FOR POLICY FUNCTIONS
OptPolInfo = struct;
OptPolInfo.policyEqNames =  {'Taylor rule';'QE rule'};
OptPolInfo.policyShockMnems =  {'etar';'etaq'};
OptPolInfo.instrumentMnems = {'ir';'qInstr'};
OptPolInfo.instrumentWeights = [0; epsilon];
OptPolInfo.objVarMnems =  {'pie';     'x';  'Dq'};
OptPolInfo.objVarWeights = [1;  pars.lambda; epsilonDelta];
OptPolInfo.beta = pars.beta;
boundInfo.instrumentMnems = {'ir'; 'qInstr';'qInstr'};
boundInfo.instrumentCoeffs = [1;      1;     -1];
boundInfo.constants =        [b;     qLB;   -qUB];
OptPolInfo.Constraints = boundInfo;

%% SPLIT STATE SPACE INTO BLOCKS (MAINLY FOR CLUSTER)
% Split state space into blocks and solve using piecewise linear method
nstilde = size(Stilde,1);
nStatesPerBlock = floor(nstilde/nBlocksForGuess);
blockStildeIndices = nan(nBlocksForGuess,2);
for iBlock = 1:nBlocksForGuess
    blockStildeIndices(iBlock,1) = (iBlock-1)*nStatesPerBlock + 1;
    if iBlock == nBlocksForGuess
        blockStildeIndices(iBlock,2) = nstilde;
    else
        blockStildeIndices(iBlock,2) = iBlock*nStatesPerBlock;
    end
end

%% DO THE LOOP
qBounds = [qLB,qUB];
x0 = zeros(nx,1);        
Shocks = struct;
Shocks.anticipated = zeros(nz,H);

polFuns0 = struct;

for iVariant = 1:2
    f0 = nan(nstilde,nEndogVars);
    Ef0 = f0;
    if iVariant == 1
        % State contingent
        iModelXiL = ModelXiL;
        iModelXiH = ModelXiH;
        condPolFunHandleQE = ...
            @(Ms,s,Ez,DEz) compute_conditional_outcomes_QE(...
            Ms,s,Ez,DEz,...
            pars.kappa,pars.beta,pars.sigma,nu,xiL,xiH,...
            pars.lambda,epsilon,epsilonDelta,b,qBounds);
    else
        % Not state contingent
        iModelXiL = ModelXiL;
        iModelXiH = ModelXiL;
        condPolFunHandleQE = ...
            @(Ms,s,Ez,DEz) compute_conditional_outcomes_QE(...
            Ms,s,Ez,DEz,...
            pars.kappa,pars.beta,pars.sigma,nu,xiL,xiL,...
            pars.lambda,epsilon,epsilonDelta,b,qBounds);
    end
    
    for iBlock = 1:nBlocksForGuess
        tic;
        disp(['Solving for block ' num2str(iBlock) ' of ' ...
            num2str(nBlocksForGuess)]);
        iBlockInds = ...
            blockStildeIndices(iBlock,1):blockStildeIndices(iBlock,2);
        iStilde = Stilde(iBlockInds,:);
        [if0,iEf0] = compute_piecewise_linear_solution_block(ModelXiL,...
            ModelXiH,OptPolInfo,iStilde);
        f0(iBlockInds,:) = if0;
        Ef0(iBlockInds,:) = iEf0;
        toc;
    end
    
    % Build a cross check f0 using conditional policy function
    OptionsOD = struct;
    OptionsOD.Algorithm.maxIter = 20000;
    OptionsOD.Algorithm.tol = 1e-06;
    %     Solve model to extract estimates of DEz from LQ problem
    BODQE = solve_LSS_model_under_optimal_discretion(ModelXiL,OptPolInfo,...
        OptionsOD);
    DEzLQxiL = BODQE(endogVarInds,lookup_model_index_numbers(xMnems,'q'))';
    BODQE = solve_LSS_model_under_optimal_discretion(ModelXiH,OptPolInfo,...
        OptionsOD);
    DEzLQxiH = BODQE(endogVarInds,lookup_model_index_numbers(xMnems,'q'))';
    [nStilde,nz] = size(f0);
    f0check = nan(nStilde,nz);
    for ijInd = 1:nStilde
        jS = Stilde(ijInd,1:2);
        jqlag = Stilde(ijInd,3);
        if jS(1)==0
            f0check(ijInd,:) = condPolFunHandleQE(jS,jqlag,...
                Ef0(ijInd,:),DEzLQxiL);
        else
            f0check(ijInd,:) = condPolFunHandleQE(jS,jqlag,...
                Ef0(ijInd,:),DEzLQxiH);
        end
    end
    
    % Pack results according to whether state contingent or not
    if iVariant == 1
        fieldName = 'stateContingent';
    else
        fieldName = 'notStateContingent';
    end
    polFuns0.(fieldName).f0 = f0;
    polFuns0.(fieldName).Ef0 = Ef0;
    polFuns0.(fieldName).condPolFunHandleQE = ...
        condPolFunHandleQE;
    polFuns0.(fieldName).f0check = f0check;
end

%% ADD ADDITIONAL PARAMETERS
pars.epsilon = epsilon;
pars.epsilonDelta = epsilonDelta;
pars.nu = nu;
pars.xiL = xiL;
pars.xiH = xiH;

%% SAVE INFORMATION REQUIRED FOR GLOBAL SOLUTION
if saveInputs
   Inputs = struct;
   Inputs.Model = Model;   
   Inputs.p = p;
   Inputs.OptPolInfo = OptPolInfo;
   Inputs.endogVarMnems = endogVarMnems;
   Inputs.Sr = Sr;
   Inputs.Svarsig = Svarsig;
   Inputs.Sq = Sq;
   Inputs.S = S;
   Inputs.Stilde = Stilde;
   Inputs.Omegavarsig = Omegavarsig;
   Inputs.Omegar = Omegar;
   Inputs.Omega = Omega;
   Inputs.condPolFunHandleQE = condPolFunHandleQE;

   Inputs.b = b;
   Inputs.qBounds = qBounds;
   Inputs.Lambda = Lambda;

   Inputs.polFuns0 = polFuns0;
   if checkConditionalOutcomes
       Inputs.f0check = f0check;
   end
   % State contingent information
   Inputs.pars = pars;
   Inputs.m = m;
   Inputs.xiL = xiL;
   Inputs.xiH = xiH;

   % Save the file
   fullSaveFileName = [saveFileFolder '\' saveFileName '.mat'];
   save_content_of_structure_to_mat_file(fullSaveFileName,Inputs);
end

%% REPORT TIME TAKEN
clockFinishTime = clock;
compute_and_display_elapsed_time(clockStartTime,clockFinishTime);