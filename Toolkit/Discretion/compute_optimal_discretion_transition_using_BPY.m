function [xtildef,rf,muf,J,diagnostics] = ...
    compute_optimal_discretion_transition_using_BPY(xtilde0,ztilde,...
    HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,beta,W,Q,S,b,J0,VxtildextildeHplus1,...
    BxtildextildeHplus1,BrxtildeHplus1,Options)
% Solves for a transition along which instrument bound constraints may bind
% Uses the BPY logic and algorithm.
%
% INPUTS:   
%   -> xtilde0: n_{\tilde{x}}*1 vector of intial conditions
%   -> ztilde: n_{\tilde{z}}*H matrix of anticipated shocks
%   -> HtildeBxtilde: HB with policy rule row and instrument(s) col removed
%   -> HtildeCxtilde: **Either**
%        (a) HC with policy rule row and instrument(s) col removed
%        (b) An n_{\tilde{x}}*n_{\tilde{x}}*H array consisting of the
%        relevant HCtildeCtilde matrix for each period t=1,...,H.
%       Note that (b) is an input when there are occasionally binding
%       constraints on non-policy variables as well as instruments.
%   -> HtildeCr: HC with policy rule row and non instrument(s) col removed
%   -> HtildeFxtilde: HF with policy rule row and instrument(s) col removed
%   -> HtildeFr: HF with policy rule row and non instrument(s) col removed
%   -> PSItildeztilde: PSI with policy rule row and shock col removed
%   -> beta: discount factor
%   -> W: target weights matrix
%   -> Q: instrument weights matrix
%   -> S: coefficients on instruments in the constraints (as in S*r>=b)
%   -> b: bounds in the inequality constraints
%   -> J0: n_{\mu}*n_{\mu}*H guess at which of the constraints is
%      binding in which period
%   -> VxtildextildeHplus1: quadratic coefficients on endogenus vars in 
%      terminal period value function
%   -> BxtildextildeHplus1: terminal regime law of motion transition for
%      non-instrument enedogneous variables
%   -> BrxtildeHplus1: terminal regime law of motion transition for 
%      instruments
%   -> Options (optional): structure containing all or any of the following
%       - tol: tolerance for convergence
%       - maxiter: maximum number before giving up
%       - reportProgress: false/true
%       - useBestOnNonConvergence: false/true
%       - issueCyclingWarning: true/false
%
% OUTPUTS:  
%   -> xtildef: n_{\tilde{x}}*H matrix of projectionss for non-instuments
%      endogenous variables
%   -> rf: n_{r}*H matrix of projectionss for instruments
%   -> muf: n_{\mu}*H matrix of lagrange multipliers
%   -> J: n_{\mu}*n_{\mu}*H indicator for which of the constraints is
%      binding in each period
%
% DETAILS: 
%   -> This function solves for a transition to a terminal OD regime under
%      the assumption that none of the bound constraints is binding in that
%      regime.
%   -> The solution algorithm is from Brendan, Paustian and Yates' paper
%      called "Optimal conventional and unconventional monetary policy in
%      the presence of collateral constraints and the zero bound".  It is
%      described in detail in Section 5 of our paper "Optimal policy with 
%      occasionally binding constraints: piecewise linear solution methods".
%   -> The algorithm also allows for anticipated shocks as described in 
%      Section 4 of the paper.
%
% NOTES: 
%   -> Error handling is kept to a minimum.
%   -> It should, in principle, be relatively straightforward to extend
%      the documentation and this code to allow for transition to a regime
%      in which one or more of the constraints is binding in perpetuity.
%
% This version: 22/12/2020
% Author(s): Matt Waldron & Richard Harrison

%% CHECK INPUTS
if nargin < 17
    error('This function requires 17 or 18 inputs')
end

%% HANDLE OPTIONAL OPTIONS INPUT
DefaultOptions = get_default_algorithm_options();
if nargin > 17
    Options = overlay_default_structure(Options,DefaultOptions);
else
    Options = DefaultOptions;
end

%% UNPACK ALGORITHM OPTIONS
maxIter = Options.maxIter;
constraintTol = Options.constraintTol;
reportProgress = Options.reportProgress;
useBestOnNonConvergence = Options.useBestOnNonConvergence;
issueCyclingWarning = Options.issueCyclingWarning;

%% COMPUTE DIMENSIONS OF PROBLEM
% It would be straightforward to do some input dimension compatibility
% checks here.
nr = size(Q,1);
nxtilde = size(W,1);
nztilde = size(PSItildeztilde,2);
nmu = size(J0,1);
H = size(ztilde,2);

%% CHECK WHETHER NON-POLICY OBCS ARE IN PLAY
if size(HtildeCxtilde,3)>1
    if size(HtildeCxtilde,3)~=H
        error(['When non-policy OBCs are applied, input "HtildeCxtilde"',...
            ' must have dimension n_{\tilde{x}}*n_{\tilde{x}}*H.']);
    end
    applyNonPolicyOBCs = true;
else
    applyNonPolicyOBCs = false;
end

%% INITIALISE SOLUTION MATRICES
Brxtilde = NaN*ones(nr,nxtilde,H);
XIrztilde = NaN*ones(nr,nztilde,H,H);
gammar = NaN*ones(nr,1,H);
Bxtildextilde = NaN*ones(nxtilde,nxtilde,H);
XIxtildeztilde = NaN*ones(nxtilde,nztilde,H,H);
gammaxtilde = NaN*ones(nxtilde,1,H);
Bmuxtilde = NaN*ones(nmu,nxtilde,H);
XImuztilde = NaN*ones(nmu,nztilde,H,H);
gammamu = NaN*ones(nmu,1,H);

%% INITIALISE CONSTANTS ASSOCIATED WITH LAGRANGE MULTIPLIERS (CONSTRAINTS)
VxtildegammaHplus1 = zeros(nxtilde,1);
gammaxtildeHplus1 = zeros(nxtilde,1);
gammarHplus1 = zeros(nr,1);

%% INITIALISE COEFFICIENTS ASSOCIATED WITH SHOCKS
% There are no shocks beyond period H so these are empty matrices.
VxtildeztildeHplus1 = zeros(nxtilde,nztilde,0);
XIxtildeztildeHplus1 = zeros(nxtilde,nztilde,0);
XIrztildeHplus1 = zeros(nr,nztilde,0);

%% INITIALISE HISTORIES
Jhist = NaN*ones(nmu,nmu,H,maxIter);
mufhist = NaN*ones(nmu,H,maxIter);
nufhist = NaN*ones(nmu,H,maxIter);
xtildefhist = NaN*ones(nxtilde,H,maxIter);
rfhist = NaN*ones(nr,H,maxIter);

%% INITIALISE ALGORITHM 
hasConverged = false;
cyclingDetected = false;
iter = 0;
Jnew = J0;

%% ALGORITHM
while ~hasConverged && iter<maxIter && ~cyclingDetected
    %% INCREASE ITERATION COUNT
    iter = iter+1;
    J = Jnew;   
    %% INITIALISE LEAD COEFFICIENTS
    BxtildextildehPlus1 = BxtildextildeHplus1;
    BrxtildehPlus1 = BrxtildeHplus1;
    XIxtildeztildehPlus1 = XIxtildeztildeHplus1;
    XIrztildehPlus1 = XIrztildeHplus1;
    gammaxtildehPlus1 = gammaxtildeHplus1;
    gammarhPlus1 = gammarHplus1;
    VxtildextildehPlus1 = VxtildextildeHplus1;
    VxtildeztildehPlus1 = VxtildeztildeHplus1;
    VxtildegammahPlus1 = VxtildegammaHplus1;
    %% RECURSE BACKWARDS
    for h = H:-1:1 
        %% COMPUTE COEFFICIENTS IN PERIOD h LAWS OF MOTION
        if applyNonPolicyOBCs
            HtildeCxtildeh = HtildeCxtilde(:,:,h);
        else
            HtildeCxtildeh = HtildeCxtilde;
        end
        [Bxtildextilde(:,:,h),Brxtilde(:,:,h),Bmuxtilde(:,:,h),...
            XIxtildeztilde(:,:,h:H,h),XIrztilde(:,:,h:H,h),...
            XImuztilde(:,:,h:H,h),gammaxtilde(:,:,h),gammar(:,:,h),...
            gammamu(:,:,h)] = compute_coefficients_in_laws_of_motion(...
            HtildeBxtilde,HtildeCxtildeh,HtildeCr,HtildeFxtilde,...
            HtildeFr,PSItildeztilde,beta,W,Q,S,b,J(:,:,h),...
            BxtildextildehPlus1,BrxtildehPlus1,XIxtildeztildehPlus1,...
            XIrztildehPlus1,gammaxtildehPlus1,gammarhPlus1,...
            VxtildextildehPlus1,VxtildeztildehPlus1,VxtildegammahPlus1,...
            nxtilde,nr,nmu,nztilde,H-h);       
        %% COMPUTE PERIOD h VALUE FUNCTION COEFFICIENTS
        % Note that these are not stored (because they are not needed to
        % compute the projection at the end) -- it would be straightforward
        % to store them and add them to the set of outputs.
        [Vxtildextilde,Vxtildeztilde,Vxtildegamma] = ...
            compute_value_function_coefficients(...
            beta,W,Q,Bxtildextilde(:,:,h),Brxtilde(:,:,h),...
            XIxtildeztilde(:,:,h:H,h),XIrztilde(:,:,h:H,h),...
            gammaxtilde(:,:,h),gammar(:,:,h),VxtildextildehPlus1,...
            VxtildeztildehPlus1,VxtildegammahPlus1,nxtilde,nztilde,H-h);
        %% UPDATE LEAD COEFICIENTS
        BxtildextildehPlus1 = Bxtildextilde(:,:,h);
        BrxtildehPlus1 = Brxtilde(:,:,h);
        XIxtildeztildehPlus1 = XIxtildeztilde(:,:,h:H,h);
        XIrztildehPlus1 = XIrztilde(:,:,h:H,h);
        gammaxtildehPlus1 = gammaxtilde(:,:,h);
        gammarhPlus1 = gammar(:,:,h);
        VxtildextildehPlus1 = Vxtildextilde;
        VxtildegammahPlus1 = Vxtildegamma;
        VxtildeztildehPlus1 = Vxtildeztilde;
    end
    %% COMPUTE PROJECTION
    [xtildef,rf,muf] = project(xtilde0,ztilde,Bxtildextilde,...
        Brxtilde,Bmuxtilde,XIxtildeztilde,XIrztilde,XImuztilde,...
        gammaxtilde,gammar,gammamu,nxtilde,nr,nmu,H);
    %% UPDATE GUESS AT BINDING CONSTRAINT INDICATOR
    [Jnew,hasConverged,nuf] = ...
        update_binding_constraint_indicator(...
        J,rf,muf,S,b,H,constraintTol,reportProgress,iter);
    %% STORE PROJECTION, CONSTRAINT CHECK AND INDICATOR FUNCTION
    Jhist(:,:,:,iter) = J;
    mufhist(:,:,iter) = muf;
    nufhist(:,:,iter) = nuf;
    xtildefhist(:,:,iter) = xtildef;
    rfhist(:,:,iter) = rf;
    %% CHECK FOR CYCLING
    % If current guess for indicator functions have been tried before and
    % were not an equilibrium, the guess update process exhibits cycles.
    if iter>1
        Jmat = reshape(Jhist(:,:,:,1:iter),H*nmu*nmu,iter);
        JtoTest = Jmat(:,1:iter-1);
        Jcurrent = Jmat(:,iter);
        Jcheck = (JtoTest==Jcurrent);
        JcheckAll = all(Jcheck);
        [~,JrepeatInds] = find(JcheckAll==true,1,'last');
        if ~isempty(JrepeatInds)
            cyclingDetected = true;
        end
    end
end

%% DEAL WITH NON-CONVERGENCE CASE
if ~hasConverged
    if ~useBestOnNonConvergence
        if cyclingDetected 
            errorReason = ['Cycling detected on iteration ' ...
                num2str(iter) '. Since useBestOnNonConvergence=false, '...
                'the search is terminated at this point.'];
        else
            errorReason = ['Maximum number of iterations exceeded without ',...
            'achieving convergence of the constraint-binding ',...
            'indicator function'];
        end
        error(errorReason)
    else
        if iter>maxIter
            warning(['Maximum number of iterations exceeded without ',...
                'achieving convergence of the constraint-binding ',...
                'indicator function ... selecting transition closest ',...
                'to constituting an equilibrium.'])
        elseif cyclingDetected && issueCyclingWarning
            warning(['Cycling detected, exiting on iteration ' ...
                num2str(iter) ' ... selecting transition closest ',...
                'to constituting an equilibrium.']);
        end
        itersWithNoConstraintViolations = ...
            find(all(all(nufhist>-constraintTol,1),2));
        if isempty(itersWithNoConstraintViolations)
            error(['There were no iterations in which none of the ',...
                'bound constraints were violated'])
        end
        [~,indOfIterWithSmallestLMviolation] = max(sum(sum(...
            (mufhist(:,:,itersWithNoConstraintViolations)<0).*...
            mufhist(:,:,itersWithNoConstraintViolations),1),2));
        bestIter = itersWithNoConstraintViolations(...
            indOfIterWithSmallestLMviolation);
        xtildef = xtildefhist(:,:,bestIter);
        rf = rfhist(:,:,bestIter);
        muf = mufhist(:,:,bestIter);
        J = Jhist(:,:,:,bestIter);     
    end
end

%% CHECK CONVERGENCE TO A STEADY STATE WITH NO BINDING CONSTRAINTS
if any(diag(J(:,:,H))==1)
    error(['Transition has not converged to a steady state in which ',...
        'no constraints are binding: constraint-binding indicator ',...
        'function has at least one unit entry in the final period of ',...
        'the projection'])
end

%% PACK DIAGNOSTICS
diagnostics = struct;
diagnostics.iter = iter;
diagnostics.cyclingDetected = cyclingDetected;
diagnostics.hasConverged = hasConverged;

end

%% FUNCTION TO COMPUTE COEFFICIENTS IN LAWS OF MOTION IN ARBITRARY PERIOD
function [Bxtildextildeh,Brxtildeh,Bmuxtildeh,XIxtildeztildeh,...
    XIrztildeh,XImuztildeh,gammaxtildeh,gammarh,gammamuh] = ...
    compute_coefficients_in_laws_of_motion(HtildeBxtilde,HtildeCxtilde,...
    HtildeCr,HtildeFxtilde,HtildeFr,PSItildeztilde,beta,W,Q,S,b,Jh,...
    BxtildextildehPlus1,BrxtildehPlus1,XIxtildeztildehPlus1,...
    XIrztildehPlus1,gammaxtildehPlus1,gammarhPlus1,VxtildextildehPlus1,...
    VxtildeztildehPlus1,VxtildegammahPlus1,nxtilde,nr,nmu,nztilde,Hztilde)
% Computes the coefficients in laws of motion for arbitrary period h.
%
% INPUTS:
%   -> HtildeBxtilde: \widetilde{H}^{B}_{\tilde{x}}
%   -> HtildeCxtilde: \widetilde{H}^{C}_{\tilde{x}}
%   -> HtildeCr: \widetilde{H}^{C}_{r}
%   -> HtildeFxtilde: \widetilde{H}^{F}_{\tilde{x}}
%   -> HtildeFr: \widetilde{H}^{F}_{r}
%   -> PSItildeztilde: \widetilde{\Psi}_{\tilde{x}\tilde{z}}
%   -> beta: discount factor
%   -> W: target weights matrix
%   -> Q: instrument weights matrix
%   -> S: coefficients in S*r>=b
%   -> b: constants in S*r>=b
%   -> Jh: \mathbb{J}_h
%   -> BxtildextildehPlus1: B_{\tilde{x}\tilde{x},h+1}
%   -> BrxtildehPlus1: B_{r\tilde{x},h+1}
%   -> XIxtildeztildehPlus1: {\Xi_{s,\tilde{x}\tilde{z},h+1}}_{s=0}^{H-h}
%   -> XIrztildehPlus1: {\Xi_{s,r\tilde{z},h+1}}_{s=0}^{H-h}
%   -> gammaxtildehPlus1: \gamma_{\tilde{x},h+1}
%   -> gammarhPlus1: \gamma_{r,h+1}
%   -> VxtildextildehPlus1: V_{\tilde{x}\tilde{x},h+1}
%   -> VxtildeztildehPlus1: {V_{s,\tilde{x}\tilde{z},h+1}}_{s=1}^{H-h}
%   -> VxtildegammahPlus1: V_{\tilde{x}\gamma,h+1}
%   -> nxtilde: number of non-instrument endogenous variables
%   -> nr: number of instruments
%   -> nmu: number of constraints
%   -> nztilde: number of shocks
%   -> Hztilde: number of periods for future (anticipated) shocks (H-h)
%
% OUTPUTS:
%   -> Bxtildextildeh: B_{\tilde{x}\tilde{x},h}
%   -> Brxtildeh: B_{r\tilde{x},h}
%   -> Bmuxtildeh: B_{\mu\tilde{x},h}
%   -> XIxtildeztildeh: {\Xi_{s,\tilde{x}\tilde{z},h}}_{s=0}^{Hztilde}
%   -> XIrztildeh: {\Xi_{s,r\tilde{z},h}}_{s=0}^{Hztilde}
%   -> XImuztildeh: {\Xi_{s,\mu\tilde{z},h}}_{s=0}^{Hztilde}
%   -> gammaxtildeh: \gamma_{\tilde{x},h}
%   -> gammarh: \gamma_{r,h}
%   -> gammamuh: \gamma_{\mu,h}

%% COMPUTE \Theta_h & \Theta_h^{-1}
Thetah = HtildeCxtilde...
    +HtildeFxtilde*BxtildextildehPlus1+HtildeFr*BrxtildehPlus1;
ThetahInv = eye(nxtilde,nxtilde)/Thetah;

%% COMPUTE \zeta_h
zetah = (ThetahInv*HtildeCr)'*(W+beta*VxtildextildehPlus1)*ThetahInv;

%% COMPUTE \Delta_{r,h}, \Delta_{r,h}^{-1}, \Delta_{x,h} & \Delta_{c,h}
Deltarh = Q+zetah*HtildeCr;
DeltarhInv = eye(nr,nr)/Deltarh;
Deltaxh = -zetah*HtildeBxtilde;
Deltach = (ThetahInv*HtildeCr)'*beta*VxtildegammahPlus1-...
    zetah*(HtildeFxtilde*gammaxtildehPlus1+HtildeFr*gammarhPlus1);

%% COMPUTE {\Delta_{z_{s},h}}_{s=0}^{H-h}
% Note that Deltaz(:,:,s) = \Delta_{z_{s-1},h} and that
% XIxtildeztildeLead(:,:,s) = \Xi_{s-1,\tilde{x}\tilde{z},h+1}, but that
% VxtildeztildeLead(:,:,s) = V_{s,\tilde{x}\tilde{z},h+1}
Deltazh = zeros(nr,nztilde,Hztilde+1);
Deltazh(:,:,1) = zetah*PSItildeztilde;
if Hztilde > 0
    for s = 1:Hztilde
        Deltazh(:,:,s+1) = ...
            (ThetahInv*HtildeCr)'*beta*VxtildeztildehPlus1(:,:,s)-...
            zetah*(HtildeFxtilde*XIxtildeztildehPlus1(:,:,s)+...
            HtildeFr*XIrztildehPlus1(:,:,s));
    end
end

%% COMPUTE \Gamma_{\mu\mu,h}, \Gamma_{rr,h}, \Gamma_{r\mu,h} & 
% \Gamma_{\mu r,h}
Gammamumuh = eye(nmu,nmu)/(eye(nmu,nmu)-Jh+Jh*S*DeltarhInv*S');
Gammarrh = DeltarhInv-DeltarhInv*S'*Gammamumuh*Jh*S*DeltarhInv;
Gammarmuh = DeltarhInv*S'*Gammamumuh;
Gammamurh = -Gammamumuh*Jh*S*DeltarhInv;

%% COMPUTE B_{r\tilde{x},h}, B_{\mu\tilde{x},h} & B_{\tilde{x}\tilde{x},h}
Brxtildeh = Gammarrh*Deltaxh;
Bmuxtildeh = Gammamurh*Deltaxh;
Bxtildextildeh = -ThetahInv*(HtildeBxtilde+HtildeCr*Brxtildeh);

%% COMPUTE {\Xi_{s,r\tilde{z},h}}_{s=0}^{Hztilde}, 
% {\Xi_{s,\mu\tilde{z},h}}_{s=0}^{Hztilde} &
% {\Xi_{s,\tilde{x}\tilde{z},h}}_{s=0}^{Hztilde}
XIrztildeh = zeros(nr,nztilde,Hztilde+1);
XImuztildeh = zeros(nmu,nztilde,Hztilde+1);
XIxtildeztildeh = zeros(nxtilde,nztilde,Hztilde+1);
for s = 0:Hztilde
    XIrztildeh(:,:,s+1) = Gammarrh*Deltazh(:,:,s+1);
    XImuztildeh(:,:,s+1) = Gammamurh*Deltazh(:,:,s+1);
    if s == 0
        XIxtildeztildeh(:,:,1) = ThetahInv*(...
            PSItildeztilde-HtildeCr*XIrztildeh(:,:,1));
    else
        XIxtildeztildeh(:,:,s+1) = -ThetahInv*(...
            HtildeFxtilde*XIxtildeztildehPlus1(:,:,s)+...
            HtildeFr*XIrztildehPlus1(:,:,s)+HtildeCr*XIrztildeh(:,:,s+1));
    end
end

%% COMPUTE \gamma_{r,h}, \gamma_{\mu,h} & \gamma_{\tilde{x},h}
gammarh = Gammarrh*Deltach+Gammarmuh*Jh*b;
gammamuh = Gammamurh*Deltach+Gammamumuh*Jh*b;
gammaxtildeh = -ThetahInv*(...
    HtildeFxtilde*gammaxtildehPlus1+HtildeFr*gammarhPlus1...
    +HtildeCr*gammarh);

end

%% FUNCTION TO COMPUTE VALUE FUNCTION COEFFICIENTS
function [Vxtildextildeh,Vxtildeztildeh,Vxtildegammah] = ...
    compute_value_function_coefficients(...
    beta,W,Q,Bxtildextildeh,Brxtildeh,XIxtildeztildeh,...
    XIrztildeh,gammaxtildeh,gammarh,VxtildextildehPlus1,...
    VxtildeztildehPlus1,VxtildegammahPlus1,nxtilde,nztilde,Hztilde)
% Computes the coefficients in value function in an arbitrary period h.
%
% INPUTS:
%   -> beta: discount factor
%   -> W: target weights matrix
%   -> Q: instrument weights matrix
%   -> Bxtildextildeh: B_{\tilde{x}\tilde{x},h}
%   -> Brxtildeh: B_{r\tilde{x},h}
%   -> XIxtildeztildeh: {\Xi_{s,\tilde{x}\tilde{z},h}}_{s=0}^{Hztilde}
%   -> XIrztildeh: {\Xi_{s,r\tilde{z},h}}_{s=0}^{Hztilde}
%   -> gammaxtildeh: \gamma_{\tilde{x},h}
%   -> gammarh: \gamma_{r,h}
%   -> VxtildextildehPlus1: V_{\tilde{x}\tilde{x},h+1}
%   -> VxtildeztildehPlus1: {V_{s,\tilde{x}\tilde{z},h+1}}_{s=1}^{H-h}
%   -> VxtildegammahPlus1: V_{\tilde{x}\gamma,h+1}
%   -> nxtilde: number of non-instrument endogenous variables
%   -> nztilde: number of shocks
%   -> Hztilde: number of periods for future (anticipated) shocks (H-h)
%
% OUTPUTS:
%   -> Vxtildextildeh: V_{\tilde{x}\tilde{x},h}
%   -> Vxtildeztildeh: {V_{s,\tilde{x}\tilde{z},h}}_{s=1}^{Hztilde+1}
%   -> Vxtildegammah: V_{\tilde{x}\gamma,h}

%% COMPUTE V_{\tilde{x}\tilde{x},h}
Vxtildextildeh = ...
    Bxtildextildeh'*(W+beta*VxtildextildehPlus1)*Bxtildextildeh+...
    Brxtildeh'*Q*Brxtildeh;

%% COMPUTE {V_{s,\tilde{x}\tilde{z},h}}_{s=1}^{Hztilde+1}
% Note the difference in indexing between Vxtildeztilde and XIxtildeztilde/
% XIrztilde.  With the former Vxtildeztilde(:,:,s) refers to 
% V_{s,\tilde{x}\tilde{z},h}.  While with the latter, e.g., 
% XIxtildeztilde(:,:,s) refers to \Xi_{s-1,\tilde{x}\tilde{z},h}.
Vxtildeztildeh = zeros(nxtilde,nztilde,Hztilde+1);
for s = 1:Hztilde+1
    Vxtildeztildeh(:,:,s) = ...
        Bxtildextildeh'*(W+beta*VxtildextildehPlus1)*XIxtildeztildeh(:,:,s)+...
        Brxtildeh'*Q*XIrztildeh(:,:,s);
    if s > 1
        Vxtildeztildeh(:,:,s) = Vxtildeztildeh(:,:,s)+...
            Bxtildextildeh'*beta*VxtildeztildehPlus1(:,:,s-1);
    end
end

%% COMPUTE V_{\tilde{x}\gamma,h}
Vxtildegammah = ...
    Bxtildextildeh'*(W+beta*VxtildextildehPlus1)*gammaxtildeh+...
    Brxtildeh'*Q*gammarh+...
    Bxtildextildeh'*beta*VxtildegammahPlus1;

end

%% FUNCTION TO PROJECT
function [xtildef,rf,muf] = project(xtilde0,ztilde,Bxtildextilde,...
    Brxtilde,Bmuxtilde,XIxtildeztilde,XIrztilde,XImuztilde,gammaxtilde,...
    gammar,gammamu,nxtilde,nr,nmu,H)
% Computes an OD projection subject to constraints and anticipated shocks.
%
% INPUTS:
%   -> xtilde0: n_{\tilde{x}}*1 vector of intial conditions
%   -> ztilde: n_{\tilde{z}}*H matrix of anticipated shocks
%   -> Bxtildextilde: {B_{\tilde{x}\tilde{x},h}}_{h=1}^{H}
%   -> Brxtilde: {B_{r\tilde{x},h}}_{h=1}^{H}
%   -> Bmuxtilde: {B_{\mu\tilde{x},h}}_{h=1}^{H}
%   -> XIxtildeztilde: {{\Xi_{s,\tilde{x}\tilde{z},h}}_{h=1}^{H}}_{s=0}^{H-h}
%   -> XIrztilde: {{\Xi_{s,r\tilde{z},h}}_{h=1}^{H}}_{s=0}^{H-h}
%   -> XImuztilde: {{\Xi_{s,\mu\tilde{z},h}}_{h=1}^{H}}_{s=0}^{H-h}
%   -> gammaxtilde: {\gamma_{\tilde{x},h}}_{h=1}^{H}
%   -> gammar: {\gamma_{r,h}}_{h=1}^{H}
%   -> gammamu: {\gamma_{\mu,h}}_{h=1}^{H}
%   -> nxtilde: number of non-instrument endogenous variables
%   -> nr: number of instruments
%   -> nmu: number of constraints
%   -> H: forecast horizon
%
% OUTPUTS:
%   -> xtildef: {\widetilde{x}_h}_{h=1}^{H}
%   -> rf: {r_h}_{h=1}^{H}
%   -> muf: {\mu_h}_{h=1}^{H}

%% INITIALISE THE OUTPUT
xtildef = zeros(nxtilde,H);
rf = zeros(nr,H);
muf = zeros(nmu,H);

%% INTIALISE LAGGED NON-INSTRUMENT ENDOGENOUS VARIABLES
xtildefLag = xtilde0;

%% COMPUTE PROJECTION
for h = 1:H
    %% COMPUTE PERIOD h ANTICIPATED SHOCK IMPACTS
    xtildeztildeImpact = zeros(nxtilde,1);
    rztildeImpact = zeros(nr,1);
    muztildeImpact = zeros(nmu,1);
    for s = 0:H-h
        xtildeztildeImpact = xtildeztildeImpact+...
            XIxtildeztilde(:,:,h+s,h)*ztilde(:,h+s);
        rztildeImpact = rztildeImpact+XIrztilde(:,:,h+s,h)*ztilde(:,h+s);
        muztildeImpact = muztildeImpact+...
            XImuztilde(:,:,h+s,h)*ztilde(:,h+s);
    end
    %% COMPUTE PERIOD h PROJECTION
    xtildef(:,h) = Bxtildextilde(:,:,h)*xtildefLag+...
        xtildeztildeImpact+gammaxtilde(:,:,h);
    rf(:,h) = Brxtilde(:,:,h)*xtildefLag+rztildeImpact+gammar(:,:,h);
    muf(:,h) = Bmuxtilde(:,:,h)*xtildefLag+muztildeImpact+gammamu(:,:,h);
    %% UPDATE LAGGED NON-INSTRUMENT ENDOGENLOUS VARIABLES
    xtildefLag = xtildef(:,h);
end

end

%% FUNCTION TO UPDATE BINDING CONSTRAINT INDICATOR
function [Jnew,hasConverged,nuf] = update_binding_constraint_indicator(...
    J,rf,muf,S,b,H,constraintTol,reportProgress,iter)
% Updates the binding constraint indicator.
%
% INPUTS:
%   -> J: n_{\mu}*n_{\mu}*H current guess for constraint binding indicator
%   -> rf: {r_h}_{h=1}^{H}
%   -> muf: {\mu_h}_{h=1}^{H}
%   -> S: coefficients in S*r>=b
%   -> b: constants in S*r>=b
%   -> H: forecast horizon
%   -> constraintTol: tolerance of allowance for constraint violations
%   -> reportProgress: true/false
%   -> iter: iteration number
%
% OUTPUTS:
%   -> Jnew: n_{\mu}*n_{\mu}*H updated guess for constraint indicator
%   -> hasConverged: false/true
%   -> nuf: nmu*H matrix, where nuf(:,h)=S*r_h-b

%% INITIALISE NEW CONSTRAINT INDICATOR
Jnew = J;

%% COMPUTE S*r-b
nuf = S*rf-b*ones(1,H);

%% CHECK CONSTRAINT AND LAGRANGE MULTIPLIER SIGN RESTRICTIONS
constraintViolationLogicals = (nuf<-constraintTol);
multiplierViolationLogicals = (muf<-constraintTol);

%% INITIALISE CONVERGENCE LOGICAL
hasConverged = false;

%% EXIT IF AN EQUILIBRIUM HAS BEEN FOUND
if ~any(any(constraintViolationLogicals)) && ...
        ~any(any(multiplierViolationLogicals))
    hasConverged = true;
    if reportProgress
        fprintf(['Equilibrium found on iteration ',num2str(iter),': ',...
            'no constraint or Lagrange multiplier sign violations\n']);
    end
    return
end

%% ADD IN A SINGLE BINDING CONSTRAINT IF ANY CONSTRAINTS ARE VIOLATED
% Note that the function has already been exited if the indicator
% \mathbb{J} constituted an equilibrium.
if any(any(constraintViolationLogicals))
    [maxViolationEachPeriod,indOfConstraintEachPeriod] = min(nuf,[],1);
    [~,periodOfLargestViolation] = min(maxViolationEachPeriod);
    indOfConstraint = indOfConstraintEachPeriod(periodOfLargestViolation);
    Jnew(indOfConstraint,indOfConstraint,periodOfLargestViolation) = 1;
    if reportProgress
        fprintf(['At least one constraint violated on iteration ',...
            num2str(iter),': Largest single violation in period ',...
            num2str(periodOfLargestViolation),' for constraint ',...
            'number ',num2str(indOfConstraint),'\n']);
    end
    return
end

%% REMOVE A SINGLE CONSTRAINT
% Note that the function has already been exited if the indicator
% \mathbb{J} constituted an equilibrium or there were any constraint
% violations above the tolerance.
if any(any(multiplierViolationLogicals))
    [maxViolationEachPeriod,indOfMultiplierEachPeriod] = min(muf,[],1);
    [~,periodOfLargestViolation] = min(maxViolationEachPeriod);
    indOfMultiplier = indOfMultiplierEachPeriod(periodOfLargestViolation);
    Jnew(indOfMultiplier,indOfMultiplier,periodOfLargestViolation) = 0;
    if reportProgress
        fprintf(['At least one Lagrange multiplier is negative on ',...
            'iteration ',num2str(iter),': Largest single violation in ',...
            'period ',num2str(periodOfLargestViolation),' for ',...
            'multiplier number ',num2str(indOfMultiplier),'\n']);
    end
    return    
end

end

%% FUNCTION CONTAINING DEFAULT OPTIONS FOR THE ALGORITHM
function DefaultOptions = get_default_algorithm_options()
% Contains the default options for the algorithm.
%
% INPUTS:   
%   -> none
%
% OUTPUTS:  
%   -> DefaultOptions: structure
%       - maxIter: maximum number of iterations
%       - constraintTol: tolerance of allowance for constraint violations
%       - reportProgress: false/true
%       - useBestOnNonConvergence: false/true

%% CREATE DEFAULT STRUCTURE
DefaultOptions.maxIter = 1000;
DefaultOptions.constraintTol = 1e-6;
DefaultOptions.reportProgress = false;
DefaultOptions.useBestOnNonConvergence = false;
DefaultOptions.issueCyclingWarning = true;

end