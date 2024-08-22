function P = compute_model_variable_covariance_matrix(...
    PHI,B,algorithmIsPrecise)
% This module computes the var-covariance matrix of LSS model variables.
% It computes the steady state variance covariance matrix of model
% variables in linear state space (LSS) models using a choice of
% algorithms.
%
% INPUTS:   
%   -> PHI: loadings on shocks in solution
%   -> B: loadings on model variables in solution
%   -> algorithmIsPrecise (optional): true if using the precise algorithm, 
%      false otherwise
%
% OUTPUTS:  
%   -> P: covariance matrix of the model variables
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> validate_LSS_model_solution_matrices
%
% DETAILS:
%   -> This function computes the variance-covariance matrix of LSS model
%      variables in steady state.
%   -> It has two algorithms for doings so. The first uses brute-force
%      matrix inversion. This method is computationally more expensive but
%      should result in a more precise matrix. 
%   -> The second approximates the answer by solving discrete Lyapunov 
%      equations using the "dlyap" routine from the MATLAB control toolbox.
%
% NOTES:
%   -> See <> for a description of the solution to MAPS LSS models.
%
% This version: 16/05/2011
% Author(s): Alex Haberis and Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of the compulsory inputs is as expected.
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
end

%% VALIDATE SOLUTION MATRICES
% Call a routine to validate the solution matrices passed in. This routine
% will throw an exception if they are not as expected.
validate_LSS_model_solution_matrices(B,PHI);

%% CHECK OPTIONAL INPUT
% If the optional third input was passed in, check that it is a logical
% scalar as expected. If not, set the algorithm choice to false (meaning
% that the approximate and faster algorithm is chosen).
if nargin > 2
    if ~islogical(algorithmIsPrecise) || ~isscalar(algorithmIsPrecise)
        errId = ['MAPS:',mfilename,':BadInput3'];
        generate_and_throw_MAPS_exception(errId);
    end
else
    algorithmIsPrecise = false;
end

%% COMPUTE COVARIANCE MATRIX OF SHOCKS
% Whichever algorithm is chosen, the covariance matrix of the shocks is
% used.
Q = PHI*PHI';

%% COMPUTE MODEL VARIABLES COVARIANCE MATRIX
% If the algorithm was chosen as "precise", (attempt to) compute the
% covariance matrix exactly. If not, use an approximation to it using a
% MATLAB control toolbox routine.
if algorithmIsPrecise
    nx = size(B,1);
    sparseB = sparse(B);
    kSparseB = kron(sparseB,sparseB);
    kSparseBtransformed = (speye(nx^2)-kSparseB)\Q(:);
    P = reshape(kSparseBtransformed,nx,nx);
else    
    P = dlyap(B,Q);
end

end