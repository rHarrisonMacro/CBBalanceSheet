%% WORKSPACE AND PATH MANAGEMENT
tidyUpAndSetPath;

%% START CLOCK FOR TIMING INFORMATION
clockStartTime = clock;

%% OPTIONS FOR WHICH VARIANT TO SOLVE 
includeStateContingentQE = true;
saveOutputs = true;
saveFileFolder = outputsFolder;

saveFileName = 'baseline'; 
initialiseFromF0check = true;

%% SET OPTIONS
Options = struct;
Options.maxit = 50000; 
Options.tol = 1e-08;
Options.updateWeight = 0.05; 
Options.reportDiagnostics = false;
% Additional options to accelerate convergence when close to solution
Options.fasterUpdateTol = 1e-08; 

%% LOAD REQUIRED INPUTS
load([inputsFolder saveFileName '.mat'],...
    'polFuns0','Stilde','Sq','Omega','Lambda');

%% SOLVE THE MODEL
polFunctions = struct;
for iVariant = 1:2
    if iVariant == 1
        if initialiseFromF0check    
            f0 = polFuns0.stateContingent.f0check;
        else
            f0 = polFuns0.stateContingent.f0;
        end
        condPolFunHandleQE = polFuns0.stateContingent.condPolFunHandleQE;
    else
        if initialiseFromF0check    
            f0 = polFuns0.notStateContingent.f0check;
        else 
            f0 = polFuns0.notStateContingent.f0;
        end
        condPolFunHandleQE = ...
            polFuns0.notStateContingent.condPolFunHandleQE;
    end
    tic;
    [f,fbar,D] = solve_QE_model_using_PFI(f0,condPolFunHandleQE,Stilde,Sq,...
        Omega,Lambda,Options);
    toc;
    %% RUN CROSS CHECK OF SOLUTION
    [nstilde,nz] = size(f);
    fCheck = nan(nstilde,nz);
    for m = 1:nstilde
        % (i) Extract latest guesses for expectations and derivatives
        ztilde = fbar(m,:);
        Dcal = D(m,:);
        % (ii) Extract exogenous and endogenous states
        Sm = Stilde(m,1:2)';
        qm = Stilde(m,3);
        % (iii) Solve for conditional outcomes
        z = condPolFunHandleQE(Sm,qm,ztilde,Dcal);
        % (iii) Load solution into latest guess
        fCheck(m,:) = z;
    end
    if iVariant ==1
        polFunctions.stateContingent.f = f;
        polFunctions.stateContingent.fbar = fbar;
        polFunctions.stateContingent.D = D;
        polFunctions.stateContingent.fCheck = fCheck;
    else
        polFunctions.notStateContingent.f = f;
        polFunctions.notStateContingent.fbar = fbar;
        polFunctions.notStateContingent.D = D;
        polFunctions.notStateContingent.fCheck = fCheck;
    end
end

%% SAVE RESULTS
if saveOutputs
    save([saveFileFolder saveFileName '.mat'],'polFunctions','Options');
end

%% REPORT TIME TAKEN
clockFinishTime = clock;
compute_and_display_elapsed_time(clockStartTime,clockFinishTime);