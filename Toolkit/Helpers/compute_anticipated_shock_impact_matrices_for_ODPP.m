function F = compute_anticipated_shock_impact_matrices_for_ODPP(...
    H,returnIntermediateHorizons,offsetHorizonForProjection,beta,Omega,...
    B,Fps,Fpol)
% Helper function for ODPPs to compute anticipated shock impact matrices.
%
% INPUTS:
%   -> H: horizon for the F matrix
%   -> returnIntermediateHorizons: true/false
%   -> offsetHorizonForProjection: true/false
%   -> beta: discount factor
%   -> Omega: loss function weights
%   -> B: state transition
%   -> Fps: \mathcal{F}^{ps} from derivation (impact via private sector)
%   -> Fpol: \mathcal{F}^{pol} from derivation (impact via policy)
%
% OUTPUTS:  
%   -> F: matrix of forward loadings with following alternative formats:
%       - nx*nx matrix representing F_{H} if returnIntermediateHorizons is 
%         false and offsetHorizonForProjection is false
%       - nx*nx matrix representing F_{H-1} if returnIntermediateHorizons
%         is false and offsetHorizonForProjection is true 
%       - nx*nx*H matrix representing {F_{1} \dots F_{H}} if
%         returnIntermediateHorizons is true and offsetHorizonForProjection
%         is false
%       - nx*nx*H matrix representing {F_{0} \dots F_{H-1}} if
%         returnIntermediateHorizons is true and offsetHorizonForProjection
%         is true
%
% DETAILS: 
%   -> This function computes the time-varying anticipated shock impact
%      matrices for optimal discretionary policy projections based on the 
%      MAPS implementation of the Dennis algorithm.
%   -> Note that the anticipated shock impact matrices cannote be derived
%      in the same way as in models solved in MAPS with non-optimal policy.
%      In particular, F_{s} ~= F^s as in the standard case (and hence the
%      output to this function is a horizon varying F matrix). This
%      reflects that the anticipated shocks affect policy optimisation via
%      their effect on future losses.
%
% NOTES: 
%   -> Error handling is minimal and is designed to guard against the
%      function handle created in 
%      "solve_LSS_model_under_optimal_discretion", which hardwires inputs
%      4-8 into the call.
%
% This version: 05/03/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 8
    error('This function requires 8 inputs');
elseif ~is_positive_real_integer(H)
    error(['1st input must be a positive real integer for the horizon ',...
        'over which to compute and return the forward loadings matrix']);
elseif ~is_logical_scalar(returnIntermediateHorizons)
    error(['2nd input must be a logical scalar (false/true) ',...
        'indicating whether to return F matrices for all horizons ',...
        'from period 1 to H or whether to return just period H']);
elseif ~is_logical_scalar(offsetHorizonForProjection)
    error(['3rd input must be a logical scalar (false/true) ',...
        'indicating whether to return F matrices for all horizons ',...
        'from period 1 to H or whether to return just period H']);    
end

%% COMPUTE F RECURSIVELY
nx = size(B,1);
FincHorizonZero = NaN(nx,nx,H+1);
SigmaIncHorizonZero = NaN(nx,nx,H+1);
FincHorizonZero(:,:,1) = eye(nx);
SigmaIncHorizonZero(:,:,1) = zeros(nx,nx);
for h = 1:H
    SigmaIncHorizonZero(:,:,h+1) = beta*B'*SigmaIncHorizonZero(:,:,h)+...
        Omega*FincHorizonZero(:,:,h);
    FincHorizonZero(:,:,h+1) = ...
        Fps*FincHorizonZero(:,:,h)+Fpol*SigmaIncHorizonZero(:,:,h+1);
end

%% EXTRACT F OUTPUT DEPENNDING ON CASE
if returnIntermediateHorizons
    if offsetHorizonForProjection
        F = FincHorizonZero(:,:,1:H);
    else
        F = FincHorizonZero(:,:,2:H+1);
    end
else
    if offsetHorizonForProjection
        F = FincHorizonZero(:,:,H);
    else
        F = FincHorizonZero(:,:,H+1);
    end
end

end