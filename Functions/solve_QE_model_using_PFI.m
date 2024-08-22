function [f,fbar,D,success] = ...
    solve_QE_model_using_PFI(fguess,condOutcomesFunHandle, ...
        Stilde,Sq,Omega,Lambda,Options)

%% 0. EXTRACT OPTIONS & INITIALISE ALGORITHM
maxit = Options.maxit;
tol = Options.tol;
updateWeight = Options.updateWeight; 
if isfield(Options,'reportDiagnostics')
    reportDiagnostics = Options.reportDiagnostics;
else
    reportDiagnostics = false;
end    
% Initialisation 
maxDiff = inf;
minMaxDiff = maxDiff;
it = 0;
% Gather and compute information
nstilde = size(Stilde,1);
nq = size(Sq,1);
ns = nstilde/nq;
qPrimeMin = Sq(1);
qPrimeMax = Sq(nq);
hq = (qPrimeMax - qPrimeMin)/(nq-1);
nz = size(fguess,2);

success = true;

%% OPTIONS TO SPEED SOLUTION CONVERGENCE
if isfield(Options,'fasterUpdateTol')
    baseUpdateWeight = updateWeight;
    adjustUpdateWeight = true;
else
    adjustUpdateWeight = false;
end

%% 1. INITIALISE GUESSES
fjMinus1 = fguess;
fj = fjMinus1;
fbarS = zeros(nstilde,nz);
fbar = zeros(nstilde,nz);
DS = zeros(nstilde,nz);
D = zeros(nstilde,nz);

%% 2. DO THE ITERATIONS
while maxDiff > tol && it<=maxit
    it = it+1;
    if it >= maxit
       success = false;
    end
    % (a) update the guess for static expectations
    for ip = 1:nq
        ipInds = ns*(ip-1)+1:ns*ip;
        fbarS(ipInds,:) = Omega*fjMinus1(ipInds,:);
    end
    % (b) Extract the vector of qhatprime values
    qPrimeVec = fjMinus1(:,3);
    % (c) Compute indices and weights of dhatVec in Sd    
    [Upsilon,Phi] = compute_1D_lin_interp_inds_and_weights(...
                        qPrimeVec,qPrimeMin,qPrimeMax,nq);
    % (d) Compute expectations
    for m = 1:nstilde
        iota1 = Upsilon(m,1);
        iota1tilde = Lambda(m,iota1);
        iota2 = Upsilon(m,2);
        iota2tilde = Lambda(m,iota2);
        fbar(m,:) = ...
            Phi(m,1)*fbarS(iota1tilde,:) + Phi(m,2)*fbarS(iota2tilde,:);
    end
    % (e) Update estimate of static derivatives
    for m = 1:nstilde
        if m <= ns
            DS(m,:) = 1/hq*(fbarS(m+ns,:) - fbarS(m,:));
        elseif m >= nstilde-ns+1
            DS(m,:) = 1/hq*(fbarS(m,:) - fbarS(m-ns,:));
        else
            DS(m,:) = 1/(2*hq)*(fbarS(m+ns,:) - fbarS(m-ns,:));
        end
    end
    % (f) Update estimates of derivatives
    for m = 1:nstilde
        iota1 = Upsilon(m,1);
        iota1tilde = Lambda(m,iota1);
        iota2 = Upsilon(m,2);
        iota2tilde = Lambda(m,iota2);
        D(m,:) = Phi(m,1)*DS(iota1tilde,:) + Phi(m,2)*DS(iota2tilde,:);
    end
    % (g) Update estimate of policy function
    for m = 1:nstilde
        % (i) Extract latest guesses for expectations and derivatives
        ztilde = fbar(m,:);
        Dcal = D(m,:);
        % (ii) Extract exogenous and endogenous states 
        Sm = Stilde(m,1:2)';
        qm = Stilde(m,3);
        % (iii) Solve for conditional outcomes
        z = condOutcomesFunHandle(Sm,qm,ztilde,Dcal);
        % (iii) Load solution into latest guess
        fj(m,:) = z;
    end
    if reportDiagnostics
        disp('Endogenous variable max abs and avg diffs:');
        disp(max(abs(fj-fjMinus1)));
        disp(mean(abs(fj-fjMinus1)));
        % Compute fraction of extrapolation and interpolation solutions
        qPrimeSol = fj(:,3);
        extrapUpFrac = sum(qPrimeSol>qPrimeMax)/nstilde;
        extrapDownFrac = sum(qPrimeSol<qPrimeMin)/nstilde;
        fprintf('Min and max QE solutions: %4.4f,%4.4f \n',...
            [min(fj(:,3)),max(fj(:,3))]);
        fprintf('Extrapolation fractions up/down: %4.4f/%4.4f \n',...
            [extrapUpFrac,extrapDownFrac]);
    end
    maxDiff = norm(fj(:)-fjMinus1(:),inf);
    minMaxDiff = min(maxDiff,minMaxDiff);
    meanDiff = mean(abs((fj(:)-fjMinus1(:))));
    % Check whether to switch to full updating
    if adjustUpdateWeight
        if maxDiff<Options.fasterUpdateTol 
            updateWeight = min(1,1.1*updateWeight);
        else
            updateWeight = baseUpdateWeight;
        end
    end
	% Update guess from previous iteration
    fjMinus1 = fjMinus1 + updateWeight*(fj-fjMinus1);
    % Print progress
    fprintf('Iteration %3.0f, norm %10.10f, mean %10.10f, min (norm) %10.10f \n',...
        [it,maxDiff,meanDiff,minMaxDiff]);
end

%% LOAD THE RESULT
f = fj;

end

