function [Bxtildextilde,Brxtilde,PHIxtildeztilde,PHIrztilde,...
    HhatCxtilde,HhatCr,Vxtildextilde,Theta,zeta,Deltar] = ...
    compute_optimal_discretion_solution_using_VFI(...
    HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,beta,W,Q,V0xtildextilde,B0xtildextilde,B0rxtilde,...
    AlgorithmOpts)
% Solves optimal discretion problem using a MAPS implementation of VFI.
% The solution is of the form:
% r_{t} = B_{r\tilde{x}}*\widetilde{x}_{t-1}...
%   +\Phi_{r\tilde{z}}*\widetilde{z}_{t}
% \widetilde{x}_{t} = B_{\tilde{x}\tilde{x}}*\widetilde{x}_{t-1}...
%   +\Phi_{\tilde{x}\tilde{z}}*\widetilde{z}_{t}
% And the fixed point for the value function coefficients is:
% V_{\tilde{x}\tilde{x}} = B_{\tilde{x}\tilde{x}}'...
%   *(W+beta*V_{\tilde{x}\tilde{x}})*B_{\tilde{x}\tilde{x}}...
%   +B_{r\tilde{x}}'*Q*B_{r\tilde{x}}
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
%   -> V0xtildextilde: initialisation/guess for value func coeffs
%   -> B0xtildextilde: initialisation for transition for endog vars
%   -> B0rxtilde: initialisation for transition for instruments
%   -> AlgorithmOpts (optional): structure with:
%       - tol: tolerance for convergence
%       - maxIter: maximum number of iterations
%
% OUTPUTS:  
%   -> Bxtildextilede: transition for endogenous variables
%   -> Brxtilde: transition for instruments
%   -> PHIxtildeztilde: shock impacts on endogenous variables
%   -> PHIrztilde: shock impacts on instruments
%   -> HhatCxtilde: coefficient on endogenous variables in targeting rules
%   -> HhatCr: coefficient on instruments in targeting rules
%   -> Vxtildextilde: coefficients on quadratic of contemporaneous 
%      endogenous variables in period-ahead losses
%   -> Theta: "solved out" coefficient on contemporaneous endogenous vars 
%      in structural equations 
%      (HtildeCx+HtildeFxtilde*Bxtildextilde+HtildeFr*Brxtilde)
%   -> zeta: coefficient common to loadings on all variables in the
%      targeting rule with the constraint substituted in
%   -> Deltar: coefficient on policy rate in rearranged targeting rule
%      with constraint substituted in
%
%
% NOTES: 
%   -> It would be straightforward to extend this function to return a
%      function handle that could be used to compute losses and loss
%      derivatives (or other intermediate outputs).
%   -> Error handling is very, very light.
%   -> It would be straightforward to amend this function to include an
%      option for damping the update.
%
% This version: 09/02/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 12
    error('This function requires 12 or 13 inputs')
end

%% GET ALGORITHM OPTIONS
DefaultAlgorithmOpts = get_default_algorithm_options();
if nargin < 13
    AlgorithmOpts = DefaultAlgorithmOpts;
else
    AlgorithmOpts = overlay_default_structure(...
        AlgorithmOpts,DefaultAlgorithmOpts);
end
tol = AlgorithmOpts.tol;
maxIter = AlgorithmOpts.maxIter;

%% COMPUTE DIMENSION OF NON-INSTRUMENT ENDOGENOUS VARIABLES
% This is for the \Theta^{-1} calculation in the algorithm.
nxtilde = size(W,1);

%% INITIALISE ALGORITHM 
hasConverged = false;
iter = 0;
Vxtildextilde = V0xtildextilde;
Bxtildextilede = B0xtildextilde;
Brxtilde = B0rxtilde;

%% ITERATE ON VALUE FUNCTION UNTIL CONVERGENCE
while ~hasConverged && iter < maxIter
    %% INCREASE ITERATION COUNT
    iter = iter+1;
    %% COMPUTE \Theta & \Theta^{-1}
    Theta = HtildeCxtilde+HtildeFxtilde*Bxtildextilede+HtildeFr*Brxtilde;
    ThetaInv = eye(nxtilde,nxtilde)/Theta;
    %% COMPUTE \zeta
    zeta = (ThetaInv*HtildeCr)'*(W+beta*Vxtildextilde)*ThetaInv;
    %% COMPUTE \Delta_r and \Delta_{\tilde{x}
    Deltar = Q+zeta*HtildeCr;
    Deltaxtilde = -zeta*HtildeBxtilde;
    %% COMPUTE A NEW B_{r\tilde{x}} and B_{\tilde{x}\tilde{x}}
    BrxtildeNew = Deltar\Deltaxtilde;
    BxtildextildeNew = -ThetaInv*(HtildeBxtilde+HtildeCr*BrxtildeNew);
    %% ITERATE ONE PERIOD BACK ON V_{\tilde{x}\tilde{x}}
    VxtildextildeNew = ...
        BxtildextildeNew'*(W+beta*Vxtildextilde)*BxtildextildeNew...
        +BrxtildeNew'*Q*BrxtildeNew;
    %% TEST FOR CONVERGENCE
    if all(all(abs(VxtildextildeNew-Vxtildextilde)<tol))
        hasConverged = true;
    else
        Vxtildextilde = VxtildextildeNew;
        Brxtilde = BrxtildeNew;
        Bxtildextilde = BxtildextildeNew;
    end
end

%% CHECK CONVERGENCE
if ~hasConverged
    error(['Maximum number of iterations exceeded without achieving ',...
        'convergence of the value function to within desired tolerance'])
end

%% COMPUTE \Delta_{\tilde{z}}
Deltaztilde = zeta*PSItildeztilde;

%% COMPUTE \Phi_{r\tilde{z}} and \Phi_{\tilde{x}\tilde{z}}
PHIrztilde = Deltar\Deltaztilde;
PHIxtildeztilde = ThetaInv*(PSItildeztilde-HtildeCr*PHIrztilde);   

%% COMPUTE TARGETING RULE COEFFICIENTS
HhatCxtilde = -(ThetaInv*HtildeCr)'*(W+beta*Vxtildextilde);
HhatCr = Q;

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