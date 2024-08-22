function z = compute_conditional_outcomes_QE(...
    s,qLag,zprime,Dzprime,kappa,beta,sigma,nu,xiL,xiH,lambda,epsilon,...
    epsilonDelta,b,qBounds)

% Computes endogenous variables conditional on state & expectations.
%
% This version incorporates state contingency

% INPUTS:
%   -> s: 1*2 vector of exogenous states (spread, rstar)
%   -> qLag: lagged QE
%   -> zprime: 1*6 vector of expectated values for inflation, output gap,
%      QE, long rate, LM on IS curve and policy rate
%   -> Dzprime: 1*6 vector of partial derivatives of expected values 
%   -> kappa: slope of Phillips curve
%   -> beta: discount factor
%   -> sigma: interest elasticity of demand
%   -> gamma: elasticity of demand with respect to QE
%   -> xi: elasticity of one-period long yield wrt change in QE
%   -> delta: SS ratio of long to short term bonds
%   -> bplusd: SS ratio of debt to GDP
%   -> omegapi: weight on inflation in the loss fn
%   -> omegax: weight on output gap in the loss fn
%   -> b: ZLB scalar value
%   -> qbounds: 1*2 vector of lower and upper bounds for QE
%
% OUTPUTS:
%   -> x: 1*6 vector of conditional outcomes for inflation, output gap,
%      policy rate, QE, change in QE and marginal continuation loss

%% UNPACK STATES, EXPECTATIONS, TARGETS AND LOSS FUNCTION WEIGHTS
varsigma = s(1);
rstar = s(2);

Epie = zprime(1);
Ex = zprime(2);
Eq = zprime(3);
% State contingent multipliers
EmuxH = zprime(6);
EmuxL = zprime(7);

DEpie = Dzprime(1);
DEx = Dzprime(2);

qLB = qBounds(1);
qUB = qBounds(2);

%% IDENTIFY WHETHER WE ARE IN THE "HIGH SPREADS" STATE OR NOT
smallNumber = 1e-06;
if varsigma>smallNumber
    creditSpreadHigh = true;
else
    creditSpreadHigh = false;
end

%% SET STATE CONTINGENT VARIABLES
if creditSpreadHigh
    xi = xiH;
else
    xi = xiL;
end
Eximux = xiH*EmuxH + xiL*EmuxL;

%% SOLUTION WHEN ZLB DOES NOT BIND
mux = 0;
pie = (1/(1+kappa^2/lambda))*(beta*Epie);
x = -kappa/lambda*pie;
q = 1/(epsilon+(1+beta)*epsilonDelta)*(...
    epsilonDelta*qLag + beta*epsilonDelta*Eq  ...
    - beta*pie*DEpie - beta*sigma*Eximux);

q = min(qUB,max(q,qLB));
ir  = (1/sigma)*(Ex-x) + Epie + (nu+xi)*q - xi*qLag + rstar - varsigma;

%% SOLUTIONS FOR ENDOGENOUS VARIABLES IF ZLB DOES BIND
if ir < b
    ir = b;
    % Solve for pie, x, q and mux by solving a system of equations
    muxLoading = DEx + sigma*DEpie + sigma*nu + sigma*xi;
    M = [0,                  1,     -sigma*(nu+xi),                     0;
         1,                  -kappa,    0,                              0;
         kappa,              lambda,    0,                              1;
         beta*DEpie,          0,     epsilon+(1+beta)*epsilonDelta, -muxLoading];
      C = [Ex - sigma*(b +xi*qLag - Epie - rstar + varsigma);
        beta*Epie;
        0; 
         beta*epsilonDelta*Eq  + epsilonDelta*qLag - beta*sigma*Eximux]; 
    Z = M\C;
    pie = Z(1);
    x = Z(2);
    q = Z(3);
    mux = Z(4);
    % Apply both upper and lower bounds to q
    if q > qUB || q < qLB
        q = min(qUB,max(q,qLB));
        ir = b;
        qtilde = (nu+xi)*q - xi*qLag;
        x = Ex - sigma*(b - qtilde - Epie - rstar + varsigma);
        pie = beta*Epie+kappa*x;
        mux = -(lambda*x + kappa*pie);
    end
end

%% SET CONDITIONAL MUX VALUES
if creditSpreadHigh
    muxH = mux;
    muxL = 0;
else
    muxL = mux;
    muxH = 0;
end

%% INCLUDE CHECK VARIABLE THAT RECORDS SPREAD
% We can inspect the vector of expectations (ie the last element of the
% expectations vector) to check that expectations of varsigma are being
% computed properly according to the Markov process matrices).
varsigCheck = varsigma;

%% VECTOR OF ENDOGENOUS VARIABLES
z = [pie x q mux ir muxH muxL varsigCheck];
end