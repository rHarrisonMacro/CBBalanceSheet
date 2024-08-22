function [Bxtildextilde,Brxtilde,PHIxtildeztilde,PHIrztilde,...
    HhatCxtilde,HhatCr,Vxtildextilde,Theta,zeta,Deltar] = ...
    compute_optimal_discretion_solution_using_Dennis_algorithm(...
    HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,beta,W,Q,B0xtildextilde,B0rxtilde,AlgorithmOpts)
% Solves optimal discretion problem using a MAPS implementation of Dennis.
% The solution is of the form:
% r_{t} = B_{r\tilde{x}}*\widetilde{x}_{t-1}...
%   +\Phi_{r\tilde{z}}*\widetilde{z}_{t}
% \widetilde{x}_{t} = B_{\tilde{x}\tilde{x}}*\widetilde{x}_{t-1}...
%   +\Phi_{\tilde{x}\tilde{z}}*\widetilde{z}_{t}
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
%   -> B0xtildextilede: initialisation/guess for transition for endog vars
%   -> B0rxtilde: initialisation/guess for transition for instruments
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
% DETAILS: 
%   -> See Appendix A of "Optimal policy with occasionally binding
%      constraints: piecewise linear solution methods" for a full 
%      derivation. The notation matches the derivation, where e.g. 
%      Bxtildextilede in the code refers to B_{\tilde{x}\tilde{x}}.
%
% NOTES: 
%   -> It would be straightforward to extend this function to return a
%      function handle that could be used to compute losses and loss
%      derivatives (or other intermediate outputs).
%   -> Error handling is very, very light.
%   -> It would be straightforward to amend this function to include an
%      option for damping the update.
%   -> It would also be relatively straightforward to replace the doubling
%      algorithm used to compute V_{\tilde{x}\tilde{x}} with a modifed
%      Lyapunov equation solver (modification necessary because of the
%      discount factor).
%
% This version: 09/02/2018
% Author(s): Matt Waldron (based on initial code by Alex Haberis)

%% CHECK INPUTS
if nargin < 11
    error('This function requires 11 or 12 inputs')
end

%% GET ALGORITHM OPTIONS
DefaultAlgorithmOpts = get_default_algorithm_options();
if nargin < 12
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
Bxtildextilde = B0xtildextilde;
Brxtilde = B0rxtilde;

%% ITERATE ON LAW OF MOTION UNTIL CONVERGENCE
while ~hasConverged && iter < maxIter
    %% INCREASE ITERATION COUNT
    iter = iter+1;
    %% COMPUTE \Theta & \Theta^{-1}
    Theta = HtildeCxtilde+HtildeFxtilde*Bxtildextilde+HtildeFr*Brxtilde;
    ThetaInv = eye(nxtilde,nxtilde)/Theta;
    %% COMPUTE V_{\tilde{x}\tilde{x}} BY FIXED POINT
    Vxtildextilde = compute_infinite_sum_using_doubling_algorithm(...
        beta*Bxtildextilde',...
        Bxtildextilde'*W*Bxtildextilde+Brxtilde'*Q*Brxtilde,...
        Bxtildextilde);
    %% COMPUTE \zeta
    zeta = (ThetaInv*HtildeCr)'*(W+beta*Vxtildextilde)*ThetaInv;
    %% COMPUTE \Delta_r and \Delta_{\tilde{x}
    Deltar = Q+zeta*HtildeCr;
    Deltaxtilde = -zeta*HtildeBxtilde;
    %% COMPUTE A NEW B_{r\tilde{x}} and B_{\tilde{x}\tilde{x}}
    BrxtildeNew = Deltar\Deltaxtilde;
    BxtildextildeNew = -ThetaInv*(HtildeBxtilde+HtildeCr*BrxtildeNew);
    %% TEST FOR CONVERGENCE
    if all(all(abs(BxtildextildeNew-Bxtildextilde)<tol)) && ...
            all(all(abs(BrxtildeNew-Brxtilde)<tol))
        hasConverged = true;
    else
        Brxtilde = BrxtildeNew;
        Bxtildextilde = BxtildextildeNew;
    end
end

%% CHECK CONVERGENCE
if ~hasConverged
    error(['Maximum number of iterations exceeded without achieving ',...
        'convergence of the law of motion to within desired tolerance'])
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
% Contains sensible defaults for options in MAPS implementation of Dennis.
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
DefaultAlgorithmOpts.maxIter = 10000;

end