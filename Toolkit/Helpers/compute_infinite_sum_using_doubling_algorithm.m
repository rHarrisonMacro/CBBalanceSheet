function V = compute_infinite_sum_using_doubling_algorithm(a1,b,a2)
% Uses a simple doubling algorithm to find a fixed point for infinite sum.
%
% INPUTS:   
%   -> a1: first power term
%   -> b: constant
%   -> a2: second power term
%
% OUTPUTS:  
%   -> V: approximation to inifinite sum
%
% DETAILS:  
%   -> Problem is: V = sum_{j=0}^{\infty} [a1^j]*b*[a2^j]
%   -> This can be written recursively as: V = b+a1*V*a2
%
% NOTES:
%   -> This is based on code by Lars Hansen and Thomas Sargent.
%   -> It would be straightforward to extend it to allow for user-defined
%      inputs for the tolerance and maximum iteration number.
%
% This version: 05/04/2017
% Author(s): Matt Waldron

%% SET MAXIMUM ITERATION NUMBER AND TOLERANCE
tol = 1e-15;
maxIts = 1000000;

%% INITIALISE
alpha1 = a1;
alpha2 = a2;
V0 = b;
it = 0;
diff = inf;

%% DOUBLING ALGORITHM
while diff > tol
    if it > maxIts
        error('Maximum number of iterations exceeded without convergence');
    end
    V1 = V0+alpha1*V0*alpha2;
    diff = max(max(abs(V1-V0)));
    V0 = V1;
    alpha1 = alpha1*alpha1;
    alpha2 = alpha2*alpha2;
    it = it+1;
end

%% CHECK OUTPUT
if ~is_finite_real_two_dimensional_numeric_matrix(V1)
    error(['Convergence has been achieved but the resulting matrix ',...
        'includes one or more NaN, inf or non-real elements']);
end

%% SET OUTPUT
V = V1;

end