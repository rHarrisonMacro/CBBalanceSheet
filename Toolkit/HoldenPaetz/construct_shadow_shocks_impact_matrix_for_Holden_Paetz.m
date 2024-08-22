function M = construct_shadow_shocks_impact_matrix_for_Holden_Paetz(...
    B,PHI,F,xmuInds,zmuInds,H)
% Constructs shadow shocks impact matrix for a Holden-Paetz inversion.
%
% INPUTS:   
%   -> B: nx*nx transition matrix in law of motion for endogenous variables
%   -> PHI: nx*nz shock impact matrtix
%   -> F: nx*nx matrix of forward loadings to compute impact of anticipated
%      shadow shocks
%   -> xmuInds: nmu*1 vector of indices for endogenous variables to which
%      each of the inequality constraints applies
%   -> zmuInds: nmu*1 corresponding vector of indices for shadow/
%      anticipated shocks to implement each of the inequality constraints
%   -> H: scalar for number of periods over which constraints might
%      be binding
%
% OUTPUTS
%   -> M: (nmu*H)*(nmu*H) matrix of impact coefficients
%
% DETAILS:
%   -> The output M is an nmu*nmu block matrix of H*H impact matrices.
%       - The blocks measure the impact of shadow shock type j (in the
%         block columns) on endogenous variable i (in the block rows), and
%         so the diagonal blocks measure the impact of the shadow shock 
%         being used to implement a constraint for a particular endogenous 
%         variable on that variable.
%       - Each block matrix measures the impact of a shadow shock 
%         anticipated in period s (in the columns) on the variables in
%         period t (in the rows).
%   -> See "invert_to_find_shadow_shocks_using_Holden_Paetz" for the sister
%      function that implements the Holden and Paetz inversion to find a
%      set of anticipated shadow shocks that respects the inequality
%      constraints.
%   -> The logic of the Holden-Paetz approach is the observation that 
%      Lagrange multipliers associated with inequality constraints are 
%      isomorphic to "shadow" shocks.  By extension, a perfect-foresight
%      expectation that an inequality constraint will by binding in future 
%      is isomorphic to an anticipated shock.  (This turns out not to be
%      true in the case of optimal discretionary policy.)
%
% NOTES:   
%   -> Note that error handling is kept to a minimum, so use with care.
%
% This version: 05/03/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 6
    error('Function requires 6 inputs');
end

%% COMPUTE NUMBER OF CONSTRAINTS
nmu = size(xmuInds,1);
if size(zmuInds,1) ~= nmu
    error(['Number of endogenous variable (4th) and shock indices ',...
        '(5th) must be identical (since these correspond to the ',...
        'indices of the variables that are subject to a constraint ',...
        'and the shadow/anticipated shocks that will be used to ',...
        'implement those constraints)']) 
end

%% INITIALISE OUTPUT
M = zeros(nmu*H,nmu*H);

%% COMPUTE ANTICIPATED SHADOW SHOCK IMPACTS BLOCK-BY-BLOCK
% In what follows, MtsAll records the effect on all endogenous variables in
% the model in period t of all anticipated shocks in the model (shadow or 
% otherwise) anticipated in period s (and known from period 1 onwards).
% Mts records the effect of the endogenous variables subject to constraints
% in period t of the shadow shocks associated with those constraints 
% binding in anticipation in period s.  M then records these effects with
% indexing consistent with the Holden and Paetz inversion.
for s = 1:H
    for t = 1:H
        if t == 1
            MtsAll = F^(s-t)*PHI;
        else
            if t <= s
                MtsAll = B*MtsAll+F^(s-t)*PHI;
            else
                MtsAll = B*MtsAll;
            end
        end
        rowInds = ((1:nmu)-1)*H+t;
        colInds = ((1:nmu)-1)*H+s;
        Mts = MtsAll(xmuInds,zmuInds);
        M(rowInds,colInds) = Mts;
    end
end

end