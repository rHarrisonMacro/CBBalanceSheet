function [xf,FOD] = project_using_OD_model(BOD,PHIOD,FODhandle,x0,Shocks)
% Projects the endogenous variables of a model solved under discretion.
%
% INPUTS:
%   -> BOD: nx*nx transtion matrix
%   -> PHIOD: nx*nz matrix of shock loadings
%   -> FODhandle: function handle to compute horizon-varying F matrices
%   -> x0: nx*1 vector of initial conditions for model variables
%   -> Shocks: struture with: 
%       - anticipated (optional): nz*H matrix of shocks
%       - unanticipated (optional): nz*H matrix of shocks
%
% OUTPUTS
%   -> xf: nx*H matrix of projected model variables
%   -> FOD: nx*nx*H matrix of horizon-varying anticipated shock multipliers
%
% DETAILS:
%   -> This function requires a function handle containing information
%   required to compute the horizon specific effects of anticipated shocks.
%
% NOTES:   
%   -> Note that error handling is kept to a minimum.
%
% This version: 05/03/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 5
    error('This function requires 5 inputs')
end

%% UNPACK SHOCKS STRUCTURE
projectionIncUnanticipatedShocks = false;
projectionIncAnticipatedShocks = false;
if isfield(Shocks,'unanticipated')
    projectionIncUnanticipatedShocks = true;
    u = Shocks.unanticipated;
    H = size(u,2);
end
if isfield(Shocks,'anticipated')
    projectionIncAnticipatedShocks = true;
    a = Shocks.anticipated;
    H = size(a,2);
    returnIntermediateHorizons = true;
    offsetHorizonForProjection = true;
    FOD = FODhandle(...
        H,returnIntermediateHorizons,offsetHorizonForProjection);
end

%% SETUP OUTPUT
nx = size(BOD,1);
xf = NaN(nx,H);

%% COMPUTE PROJECTION
if ~projectionIncAnticipatedShocks && projectionIncUnanticipatedShocks
    uImpact = zeros(nx,1);
    for s = 1:H
        uImpact = BOD*uImpact+PHIOD*u(:,s);
        xf(:,s) = BOD^s*x0+uImpact;
    end
elseif projectionIncAnticipatedShocks && ~projectionIncUnanticipatedShocks
    R = compute_anticipated_shocks_impact(PHIOD,FOD,a,H,nx);
    aImpact = zeros(nx,1);
    for s = 1:H
        aImpact = BOD*aImpact+R(:,s);
        xf(:,s) = BOD^s*x0+aImpact;
    end
else
    R = compute_anticipated_shocks_impact(PHIOD,FOD,a,H,nx);
    aImpact = zeros(nx,1);
    uImpact = zeros(nx,1);
    for s = 1:H
        aImpact = BOD*aImpact+R(:,s);
        uImpact = BOD*uImpact+PHIOD*u(:,s);
        xf(:,s) = BOD^s*x0+aImpact+uImpact;
    end  
end

end

%% FUNCTION TO COMPUTE ANTICIPATED SHOCKS IMPACT
function R = compute_anticipated_shocks_impact(PHI,FOD,a,H,nx)
% This helper computes the impact of the anticipated shocks. 
%
% INPUTS:   
%   -> PHI: shock loadings
%   -> FOD: nx*nx*H matrix of horizon-varying forward loadings
%   -> a: nz*H matrix of anticipated shocks over the projection
%   -> H: horizon for the projection
%   -> nx: dimension of the endogenous variables
%
% OUTPUTS:  
%   -> R: nx*H matrix of anticipated shock impacts

%% INITIALISE THE OUTPUT
R = zeros(nx,H);

%% COMPUTE R PERIOD BY PERIOD
% Note that unlike in the standard case of model solved under non-optimal 
% policy (where F_s = F^{s}), it is necessary to compute the impacts in a
% double loop.
for h = 1:H
    for s = h:H
        R(:,h) = R(:,h)+FOD(:,:,s-h+1)*PHI*a(:,s);
    end
end

end