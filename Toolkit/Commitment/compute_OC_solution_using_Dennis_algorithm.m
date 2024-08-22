function [By,PHIy,Fy,HBy,HCy,HFy,PSIy,OCinfo,indicesInfo] = ...
    compute_OC_solution_using_Dennis_algorithm(...
    HtildeBxtilde,HtildeCxtilde,HtildeCr,HtildeFxtilde,HtildeFr,...
    PSItildeztilde,beta,OCinfo)
% This function computes the rational expectations solution of a model
% under optimal commitment policy, using the approach of Dennis (2007).
%
% INPUTS:   
%   -> HtildeBxtilde: lagged endog var loadings in partitioned model eqs
%   -> HtildeCxtilde: contemp. endog var loadings in partitioned model eqs
%   -> HtildeCr: contemp. instrument loadings in partitioned model eqs
%   -> HtildeFxtilde: lead endog var loadings in partitioned model eqs
%   -> HtildeFr: lead instrument loadings in partitioned model eqs
%   -> PSItildeztilde: non-policy shock loadings in partitioned model eqs
%   -> beta: discount factor
%   -> OCinfo: Structure containing information about policy problem:
%   	- W: weights on endogenous variables in loss function
%       - Q: weights on instruments in loss function
%       - instrumentIndices: indices of policy instruments
%       - Constraints (optional): structure containing information
%   about bounds on the policy instruments, with the following fields:
%           - instrumentMnems: single string or nr*1 column cell string array
%           - instrumentCoeffs: numeric column vector of coeffs on instruments
%           - constants: numeric column vector of constraint constants
%
% OUTPUTS:  
%   -> By: transition for endogenous variables
%   -> HBy: coeeficients on lags in structural form
%   -> HCy: coefficients on contemporaneous variables in structural form
%   -> HFy: coefficients on leads in structural form
%   -> PSIy: coefficients on shocks in structural form
%   -> OCinfo: Updated structure containing information about OC problem
%   -> indicesInfo: Structure containing information about the indices of
%   the endogenous variables, instruments and lagrange multipliers
%
% DETAILS: 
%   -> The code accepts a model in the following structural form:
%
% HtildeFxtilde*xtilde{t+1} + HtildeCxtilde*xtilde{t} +
% HtildeBxtilde*xtilde{t-1} + HtildeFr*r{t+1} + HtildeCr*r{t} =
% PSItildeztilde*ztilde{t}
% where xtilde are the non-policy endogenous variables, rtilde are the
% policy instrument(s) and ztilde are the non-policy shocks.
%
% The model computes the rational expectations solution for a vector y,
% given by:
%
%     [  x   ]
% y = |      |
%     [lambda]
%
% where x is the vector of (policy and non-policy) endogenous variables,
% ordered in the same way as the original MAPS model and lambda is the
% vector of Lagrange multipliers on the non-policy variables.
%
% The model under optimal policy has the following structural form:
%
% HFy*y{t+1} + HCy*y{t} + HBy*y{t-1} = PSIy*zhat{t}
%
% where zhat{t} a vector stacking ztilde and a set of dummy shocks used to
% impose instrument constraints if these are specified as inputs.
%
% The rational expectations solution is given by:
% 
% y{t} = By*y{t-1} + \sum_{i=0}^{\infty} Fy^{i}*PHIy*zhat{t+i}
%
%
% This version: 25/07/2019
% Author(s): Richard Harrison

%% DEAL WITH DEFAULT OPTIONS
if isfield(OCinfo,'Constraints')
    Constraints = OCinfo.Constraints;
    instrumentBoundsApply = true;
    S = create_instrument_bound_constraint_matrices(...
        OCinfo.instrumentMnems,...
        Constraints);
    nShadowShocks = size(S,1);
else
    instrumentBoundsApply = false;
end


%% UNPACK OCINFO STRUCTURE
W = OCinfo.W;
Q = OCinfo.Q;
instrumentInds = OCinfo.instrumentInds;

%% GATHER SIZE OF SYSTEM
nxtilde = size(W,1);
nr = size(Q,1);
nztilde = size(PSItildeztilde,2);
nlambda = nxtilde;
ny = nr + nxtilde + nlambda;

%% USE THIS TO COMPUTE INDICES OF THE DIFFERENT VARIABLE TYPES
nx = nxtilde + nr;
rInds = instrumentInds;
xtildeInds = 1:nx;
xtildeInds(rInds) = [];
lambdaInds = nx+1:ny;

%% BUILD SYSTEM OF FIRST ORDER CONDITIONS IN STRUCTURAL FORM
% Initialise matrices
HFy = zeros(ny,ny);
HCy = zeros(ny,ny);
HBy = zeros(ny,ny);
PSIyztilde = zeros(ny,nztilde);
% Load first order conditions for instrument:
HCy(1:nr,rInds) = Q;
HCy(1:nr,lambdaInds) = -HtildeCr';
HBy(1:nr,lambdaInds) = -1/beta*HtildeFr';
% Load first order conditions for non-policy variables:
HFy(nr+1:nx,lambdaInds) = -beta*HtildeBxtilde';
HCy(nr+1:nx,xtildeInds) = W;
HCy(nr+1:nx,lambdaInds) = -HtildeCxtilde';
HBy(nr+1:nx,lambdaInds) = -1/beta*HtildeFxtilde';
% Load first order conditions with respect to multipliers:
HFy(nx+1:ny,xtildeInds) = HtildeFxtilde;
HFy(nx+1:ny,rInds) = HtildeFr;
HCy(nx+1:ny,xtildeInds) = HtildeCxtilde;
HCy(nx+1:ny,rInds) = HtildeCr;
HBy(nx+1:ny,xtildeInds) = HtildeBxtilde;
PSIyztilde(nx+1:ny,:) = PSItildeztilde;

%% ADD "SHADOW SHOCKS" TO IMPLEMENT INSTRUMENT BOUNDS IF NECESSARY
if instrumentBoundsApply
    PSIy = zeros(ny,nztilde+nShadowShocks);
    PSIy(:,1:nztilde) = PSIyztilde;
    PSIy(1:nr,nztilde+1:nztilde+nShadowShocks) = S';
    % Augment OCinfo with bounds information
    if nShadowShocks == 1
        if iscell(Constraints.instrumentMnems)
            instrumentMnem = Constraints.instrumentMnems{:};
        else
            instrumentMnem = Constraints.instrumentMnems;
        end
        shadowShockMnems = ['shadowShock_1_' instrumentMnem];
    else
        shadowShockMnems = cell(nShadowShocks,1);
        for iShadowShock = 1:nShadowShocks
            shadowShockMnems{iShadowShock} = ['shadowShock_' ...
                num2str(iShadowShock) '_' ...
                Constraints.instrumentMnems{iShadowShock}];
        end
    end
    OCinfo.Constraints.shadowShockMnems = shadowShockMnems;
else
    PSIy = PSIyztilde;
    nShadowShocks = 0;
end

%% SOLVE THE MODEL
[By,PHIy,Fy] = solve_LSS_model_RE_dynamics(HBy,HCy,HFy,PSIy);

%% PACK INFORMATION ABOUT MODEL INDICES
indicesInfo = struct;
indicesInfo.xtildeInds = xtildeInds;
indicesInfo.rInds = rInds;
indicesInfo.lambdaInds = lambdaInds;
indicesInfo.nShadowShocks = nShadowShocks;
indicesInfo.nztilde = nztilde;
indicesInfo.nx = nx;
indicesInfo.ny = ny;
indicesInfo.nxtilde = nxtilde;
indicesInfo.nr = nr;
indicesInfo.nlambda = nlambda;

end