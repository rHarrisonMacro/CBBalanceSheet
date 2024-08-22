function [indices,weights,rawIndices] = ...
    compute_1D_lin_interp_inds_and_weights(vectorToTest,minVal,...
    maxVal,nGridPoints)
%compute_one_dim_lin_interp_inds_and_weights This function computes indices
%and weights required to compute one-dimensional linear interpolation/
%extrapolation for an evenly spaced grid.
%   
% INPUTS
% -> vectorToTest: nx*1 vector of values for which weights and indices
% should be computed
% -> minVal: real scalar lowest value on grid
% -> maxVal: real scalar highest value on grid
% -> nGridPoints: scalar number of (evenly spaced!) grid points
%
% OUTPUTS
% -> indices: nx*2 matrix of indices of grid points to use for linear
% interpolation
% -> weights: nx*2 matrix of weights to apply to grid points

% ISSUES:
% Extrapolation does not seem to generate the correct weights. Downward
% extrapolation generates weights that when multiplied by the relevant grid
% values generate the **incremental** extrapolation values.  So if the
% lowest two grid values are 1 and 2, extrapolation weights designed to
% deliver 0.5 generate -0.5 when applied to the grid point values (1 and
% 2).

%% GATHER KEY INFORMATION AND INITIALISE OUTPUTS
d = (maxVal-minVal)/(nGridPoints-1);
nVals = size(vectorToTest,1);
indices = nan(nVals,2);
weights = nan(nVals,2);

%% IDENTIFY INTERPOLATION AND EXTRAPOLATION
rawIndices = vectorToTest/d - minVal/d;
indsToExtrapolateUp = (rawIndices>nGridPoints-2);
indsToExtrapolateDown = (rawIndices<0);
indsToInterpolate = ~(indsToExtrapolateUp|indsToExtrapolateDown);

%% INTERPOLATION
interpFloor = floor(rawIndices(indsToInterpolate));
indices(indsToInterpolate,1) = interpFloor +1;
indices(indsToInterpolate,2) = interpFloor + 2;
weights(indsToInterpolate,2) = rawIndices(indsToInterpolate) - ...
    interpFloor;
weights(indsToInterpolate,1) = 1-weights(indsToInterpolate,2);

%% EXTRAPOLATION DOWNWARDS
indices(indsToExtrapolateDown,1) = 1;
indices(indsToExtrapolateDown,2) = 2;
weights(indsToExtrapolateDown,2) = rawIndices(indsToExtrapolateDown);
weights(indsToExtrapolateDown,1) = 1 - weights(indsToExtrapolateDown,2);

%% EXTRAPOLATION UPWARDS
indices(indsToExtrapolateUp,1) = nGridPoints-1;
indices(indsToExtrapolateUp,2) = nGridPoints;
weights(indsToExtrapolateUp,1) = (maxVal-minVal)/d - rawIndices(indsToExtrapolateUp);
weights(indsToExtrapolateUp,2) = 1 - weights(indsToExtrapolateUp,1);

end

