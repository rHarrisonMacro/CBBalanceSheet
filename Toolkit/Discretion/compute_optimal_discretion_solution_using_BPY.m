function [Bxtildextilde,Brxtilde,Bmuxtilde,PHIxtildeztilde,PHIrztilde,...
    PHImuztilde,gammaxtilde,gammar,gammamu,Vxtildextilde,Vxtildegamma,...
    Theta,zeta,Deltar] = ...
    compute_optimal_discretion_solution_using_BPY(...
    HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,beta,W,Q,S,b,J,V0xtildextilde,V0xtildegamma,...
    B0xtildextilde,B0rxtilde,gamma0xtilde,gamma0r,AlgorithmOpts)
% Solves for optimal discretionary steady state with bound constraints.
% The fixed points for the value function coefficients are:
% V_{\tilde{x}\tilde{x}} = B_{\tilde{x}\tilde{x}}'...
%   *(W+beta*V_{\tilde{x}\tilde{x}})*B_{\tilde{x}\tilde{x}}...
%   +B_{r\tilde{x}}'*Q*B_{r\tilde{x}}
% V_{\tilde{x}\gamma} = B_{\tilde{x}\tilde{x}}'...
%   *(W+beta*V_{\tilde{x}\tilde{x}})*\gamma_{\tilde{x}}...
%   +B_{r\tilde{x}}'*Q*\gamma_{r}...
%   +B_{\tilde{x}\tilde{x}}'*V_{\tilde{x}\gamma}
% The coefficients in the laws of motion consistent with that steady state
% correspond to the following:
% r_{t} = B_{r\tilde{x}}*\widetilde{x}_{t-1}...
%   +\Phi_{r\tilde{z}}*\widetilde{z}_{t}+\gamma_{r}
% \widetilde{x}_{t} = B_{\tilde{x}\tilde{x}}*\widetilde{x}_{t-1}...
%   +\Phi_{\tilde{x}\tilde{z}}*\widetilde{z}_{t}+\gamma_{\tilde{x}}
% \mu_{t} = B_{\mu\tilde{x}}*\widetilde{x}_{t-1}...
%   +\Phi_{\mu\tilde{z}}*\widetilde{z}_{t}+\gamma_{\mu}
%
% INPUTS:   
%   -> HtildeBxtilde: lagged endog var loadings in partitioned model eqs
%   -> HtildeCxtilde: contemp. endog var loadings in partitioned model eqs
%   -> HtildeCr: contemp. instrument loadings in partitioned model eqs
%   -> HtildeFxtilde: lead endog var loadings in partitioned model eqs
%   -> HtildeFr: lead instrument loadings in partitioned model eqs
%   -> PSItildeztilde: non-policy shock loadings in partitioned model eqs
%   -> beta: discount factor
%   -> W: weights on endogenous variables in loss function
%   -> Q: weights on instruments in loss function
%   -> S: coefficients on instruments in inequality constraints
%   -> b: constants in inequality constraint
%   -> J: indicator for which of the constraints is binding in steady state
%   -> V0xtildextilede: initialisation/guess for value func quad coeffs
%   -> V0xtildegamma: initialisation/guess for value func constant coeffs
%   -> B0xtildextilde: initialisation for transition for endog vars
%   -> B0rxtilde: initialisation for transition for instruments
%   -> gamma0xtilde: initialisation for constant for endog vars
%   -> gamma0r: initialisation for constant for instruments
%   -> AlgorithmOpts (optional): structure with:
%       - tol: tolerance for convergence
%       - maxIter: maximum number of iterations
%
% OUTPUTS:  
%   -> Bxtildextilde: transition for endogenous variables
%   -> Brxtilde: transition for instruments
%   -> Bmuxtilde: transition for lagrange multipliers
%   -> PHIxtildeztilde: shock impacts on endogenous variables
%   -> PHIrztilde: shock impacts on instruments
%   -> PHImuztilde: shock impacts on lagrange multipliers
%   -> gammaxtilde: constants on endogenous variables
%   -> gammar: constants on instruments
%   -> gammamu: constants on Lagrange multipliers
%   -> Vxtildextilde: coefficients on quadratic of contemporaneous 
%      endogenous variables in period-ahead losses
%   -> Vxtildegamma: coefficients on constants of contemporaneous 
%      endogenous variables in period-ahead losses
%   -> Theta: "solved out" coefficient on contemporaneous endogenous vars 
%      in structural equations 
%      (HtildeCx+HtildeFxtilde*Bxtildextilde+HtildeFr*Brxtilde)
%   -> zeta: coefficient common to loadings on all variables in the
%      targeting rule with the constraint substituted in
%   -> Deltar: coefficient on policy rate in rearranged targeting rule
%      with constraint substituted in
%
% DETAILS: 
%   -> See Section 5 of "Optimal policy with occasionally binding
%      constraints: piecewise linear solution methods" for a derivation. 
%      The notation matches the derivation, where e.g. Bxtildextilede in 
%      the code refers to B_{\tilde{x}\tilde{x}}.
%   -> If the indicator J (\mathbb{J} in the derivation) is a null matrix
%      (of dimension n_{\mu}  by n_{\mu}), then this function should
%      produce exactly the same result as
%      "compute_optimal_discretion_solution_using_VFI".
%
% NOTES: 
%   -> It would be straightforward to extend this function to return the 
%      constant term in the loss function and/or any intermediate
%      calculations.
%   -> Error handling is very, very light.
%   -> It would be straightforward to amend this function to include an
%      option for damping the update.
%
% This version: 26/02/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 18
    error('This function requires 18 or 19 inputs')
end

%% GET ALGORITHM OPTIONS
DefaultAlgorithmOpts = get_default_algorithm_options();
if nargin < 19
    AlgorithmOpts = DefaultAlgorithmOpts;
else
    AlgorithmOpts = overlay_default_structure(...
        AlgorithmOpts,DefaultAlgorithmOpts);
end
tol = AlgorithmOpts.tol;
maxIter = AlgorithmOpts.maxIter;

%% COMPUTE DIMENSION OF NON-INSTRUMENT ENDOGENOUS VARIABLES & CONSTRAINTS
% This is for the \Theta^{-1}, \Delta_r^{-1} and \Sigma_{\mu\mu} 
% calculations in the algorithm.
nxtilde = size(W,1);
nr = size(Q,1);
nmu = size(b,1);

%% INITIALISE ALGORITHM 
hasConverged = false;
iter = 0;
Vxtildextilde = V0xtildextilde;
Vxtildegamma = V0xtildegamma;
Bxtildextilde = B0xtildextilde;
Brxtilde = B0rxtilde;
gammaxtilde = gamma0xtilde;
gammar = gamma0r;

%% ITERATE ON VALUE FUNCTION UNTIL CONVERGENCE
while ~hasConverged && iter < maxIter
    %% INCREASE ITERATION COUNT
    iter = iter+1;
    %% COMPUTE \Theta & \Theta^{-1}
    Theta = HtildeCxtilde+HtildeFxtilde*Bxtildextilde+HtildeFr*Brxtilde;
    ThetaInv = eye(nxtilde,nxtilde)/Theta;
    %% COMPUTE \zeta
    zeta = (ThetaInv*HtildeCr)'*(W+beta*Vxtildextilde)*ThetaInv;
    %% COMPUTE \Delta_r, \Delta_r^{-1}, \Delta_{\tilde{x} & \Delta_c
    Deltar = Q+zeta*HtildeCr;
    DeltarInv = eye(nr,nr)/Deltar;
    Deltaxtilde = -zeta*HtildeBxtilde;
    Deltac = (ThetaInv*HtildeCr)'*beta*Vxtildegamma...
        -zeta*(HtildeFxtilde*gammaxtilde+HtildeFr*gammar);
    %% COMPUTE \Gamma_{\mu\mu}, \Gamma_{rr} & \Gamma_{r\mu}
    Gammamumu = eye(nmu,nmu)/(eye(nmu,nmu)-J+J*S*DeltarInv*S');
    Gammarr = DeltarInv-DeltarInv*S'*Gammamumu*J*S*DeltarInv;
    Gammarmu = DeltarInv*S'*Gammamumu;
    %% COMPUTE B_{r\tilde{x}}, B_{\tilde{x}\tilde{x}}, \gamma_{r} & \gamma_{\tilde{x}}
    BrxtildeNew = Gammarr*Deltaxtilde;
    BxtildextildeNew = -ThetaInv*(HtildeBxtilde+HtildeCr*BrxtildeNew);
    gammarNew = Gammarr*Deltac+Gammarmu*J*b;
    gammaxtildeNew = -ThetaInv*...
        (HtildeFxtilde*gammaxtilde+HtildeFr*gammar+HtildeCr*gammarNew);
    %% COMPUTE NEW V_{\tilde{x}\tilde{x}} & V_{\tilde{x}\gamma}
    VxtildextildeNew = ...
        BxtildextildeNew'*(W+beta*Vxtildextilde)*BxtildextildeNew...
        +BrxtildeNew'*Q*BrxtildeNew;
    VxtildegammaNew = ...
        BxtildextildeNew'*(W+beta*Vxtildextilde)*gammaxtildeNew...
        +BrxtildeNew'*Q*gammarNew+BxtildextildeNew'*Vxtildegamma;
    %% TEST FOR CONVERGENCE
    if all(all(abs(VxtildextildeNew-Vxtildextilde)<tol)) && ...
           all(all(abs(VxtildegammaNew-Vxtildegamma)<tol)) 
        hasConverged = true;
    else
        Vxtildextilde = VxtildextildeNew;
        Vxtildegamma = VxtildegammaNew;
        Brxtilde = BrxtildeNew;
        Bxtildextilde = BxtildextildeNew;
        gammar = gammarNew;
        gammaxtilde = gammaxtildeNew;
    end
end

%% CHECK CONVERGENCE
if ~hasConverged
    error(['Maximum number of iterations exceeded without achieving ',...
        'convergence of the value function to within desired tolerance'])
end

%% COMPUTE \Delta_{\tilde{z}}
Deltaztilde = zeta*PSItildeztilde;

%% COMPUTE \Gamma_{\mu r}
Gammamur = -Gammamumu*J*S*DeltarInv;

%% COMPUTE B_{\mu\tilde{x}} & \gamma_{\mu}
Bmuxtilde = Gammamur*Deltaxtilde;
gammamu = Gammamur*Deltac+Gammamumu*J*b;

%% COMPUTE \Phi_{r\tilde{z}}, \Phi_{\tilde{x}\tilde{z}} & \Phi_{\mu\tilde{z}}
PHIrztilde = Gammarr*Deltaztilde;
PHIxtildeztilde = ThetaInv*(PSItildeztilde-HtildeCr*PHIrztilde);   
PHImuztilde = Gammamur*Deltaztilde;

end

%% FUNCTION CONTAINING DEFAULT ALGORITHM OPTIONS
function DefaultAlgorithmOpts = get_default_algorithm_options()
% Contains sensible defaults (hopefully).
%
% INPUTS:   
%   -> none
%
% OUTPUTS:  
%   -> DefaultAlgorithmOpts: structure
%       - tol: tolerance for convergence
%       - maxIter: maximum number of iterations

%% CREATE DEFAULT OPTIONS STRUCTURE
DefaultAlgorithmOpts.tol = 1e-12;
DefaultAlgorithmOpts.maxIter = 100000;

end