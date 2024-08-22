function [B,PHI,F] = solve_LSS_model_RE_dynamics(HB,HC,HF,PSI)
% This module solves LSS models under rational expectations.
% It contains two algorithms for a linear state space (LSS) model''s 
% solution: one for backward-looking models and one for forward-looking 
% models. The algorithm used for the latter is the FED's AIM algorithm.
%
% INPUTS:   
%   -> HB: nx*nx matrix of loadings on lagged model variables
%   -> HC: nx*nx matrix of loadings on contemporaneous model variables
%   -> HF: nx*nx matrix of loadings on future expected model variables
%   -> PSI: nx*nz matrix of loadings on shocks
%
% OUTPUTS:  
%   -> B: nx*nx matrix of loadings on lagged model variables in solution
%   -> PHI: nx*nz matrix of loadings on shocks in solution
%   -> F: nx*nz matrix of forward loadings on model variables (for
%      anticipated shocks)
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> validate_LSS_model_equation_matrices
%   -> SPAmalg (to be updgraded)
%
% DETAILS:  
%   -> This module solves LSS models for their rational expectations (RE) 
%      dynamics. 
%   -> If the input model is forward-looking, it uses the AIM algorithm to
%      solve for the dynamics.
%   -> If it is backward-looking is inverts the contemporaneous loadings
%      directly.
%
% NOTES:
%   -> See <> for a description of LSS model's and their solution.
%   -> Note that the AIM algorithm has been taken directly from the
%      internet download. MAPS will be upgraded in the future to include a
%      re-written algorithm with more comprehensive error handling.
%
% This version: 16/05/2011
% Author(s): Alex Haberis and Matt Waldron

%% CHECK INPUTS
% Check that the correct number of inputs was passed in.
if nargin < 4
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
end

%% VALIDATE STRUCTURAL MATRICES
% Call a validation routine to check that the matrices are of the correct
% shape and are consistent with each other.
validate_LSS_model_equation_matrices(HB,HC,HF,PSI);

%% CHECK MODEL IS VALID
% Throw an exception of the model contains no contemporaneous terms.
if ~any(any(HC))
    errId = ['MAPS:',mfilename,':NoContemporaneousVars'];
    generate_and_throw_MAPS_exception(errId);
end

%% COMPUTE DIMENSION OF MODEL VARIABLES
% Compute the dimension of the model variables for use below.
nx = size(HC,1);

%% BACKWARD-LOOKING MODEL SOLUTION
% Solve backward-looking model by inverting the HC matrix directly. Set the
% forward loadings on shocks to zeros (because anticipated shocks are not 
% defined in backward-looking models). 
if ~any(any(HF))
    B = -HC\HB;
    PHI = HC\PSI;
    F = zeros(nx,nx);
end

%% FORWARD-LOOKING MODEL SOLUTION
% Solev forward-looking model using the AIM algorithm. Throw an exception
% if the algorithm fails. Otherwise, compute the shock loadings and forward
% loadings from the AIM output and structural matrices.
if any(any(HF))
    [B,~,~,~,~,~,aimcode] = SPAmalg([HB HC HF],nx,1,1,1e-10,1);
    if aimcode ~= 1
        errId = ['MAPS:',mfilename,':AIMunstableSolution'];
        generate_and_throw_MAPS_exception(errId);
    end
    phiMat = (HC+HF*B)\eye(nx);
    PHI = phiMat*PSI;
    F = -phiMat*HF;
end

end