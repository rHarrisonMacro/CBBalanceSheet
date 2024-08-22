function [mustar,areConstraintsBinding] = ...
    invert_to_find_shadow_shocks_using_Holden_Paetz(xf,bounds,M,Options)
% Computes shadow shocks in a Holden-Paetz inversion.
% The marginal set of anticipated shocks output can be combined with a
% baseline set of shocks to produce a projection in which one or more of
% the endogenous variables is subject to an inequality constraint.
%
% INPUTS:   
%   -> xf: nmu*T matrix containing an unconstrained baseline projection for
%      endogenous variables subject to the constraints
%   -> bounds: nmu*1 vector of bounds for the inequality constraints
%   -> M: (nmu*T)*(nmu*T) matrix of impact coefficients on anticipated
%      shadow shocks
%   -> Options (optional): structure of options for the optimisation
%       - a0 (optional): initial condition for the search (must be 
%         specified as an nmu*H matrix (rows indexing shadow shock, 
%         cols periods)
%       - fvalTol (optional): tolerance for accepting the Linear 
%         Complementarity Solution
%       - Quadprog (optional): structure of options for the QP algorithm
%         (see sub-function below for info on valid content)
%
% OUTPUTS
%   -> mustar: nmu*T matrix of anticipated shadow shocks necessary to
%      impose the constraints
%   -> areConstraintsBinding: true/false
%
% DETAILS:
%   -> This function implements a Holden and Paetz inversion to find a 
%      marginal set of anticipated shadow shocks that ensures a set of
%      inequality constraints is respected.
%   -> The approach is explained in their 2012 paper entitled "Efficient
%      simulation of DSGE models with inequality constraints".  It uses
%      quadratic programming to find a (there may exist multiple 
%      equilibria) set of shocks that solve a complementary slackness
%      condition (that summarises the constraint) and respect the
%      constraint that the shocks must be non-strictly positive.
%   -> The M matrix input records the impact of all the anticipated 
%      shadow shocks (ordered in a particular way).  See
%      "construct_shadow_shocks_impact_matrix_for_Holden_Paetz" for
%      details.
%
%
% This version: 22/03/2018
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 3
    error('Function requires 3 or 4 inputs');
end

%% COMPUTE DIMENSIONALITY OF PROBLEM
[nmu,T] = size(xf);

%% COMPUTE DISTANCE OF BASE PROJECTION FROM BOUNDS FOR QUADPROG
% Note that this is transposed because Holden and Paetz stacks variables
% and shocks into columns with time indexed first and then variables after.
xfminusBounds = xf'-ones(T,1)*bounds';
v = xfminusBounds(:);

%% CHECK CONSTRAINTS ARE VIOLATED
% If no constraints are violated, then the shadow shocks are zero (by 
% definition) and we can ext the function. 
if all(v>=0)
    areConstraintsBinding = false;
    mustar = zeros(nmu,T);
    return
else
    areConstraintsBinding = true;
end

%% HANDLE OPTIONAL OPTIONS INPUT
DefaultOpts = get_default_options(nmu,T);
mu0 = DefaultOpts.mu0;
fvalTol = DefaultOpts.fvalTol;
QuadProgOpts = DefaultOpts.QuadProg;
if nargin > 3
    check_field_names_in_structure(...
        Options,'options','',{'a0';'fvalTol';'Quadprog'});
    if isfield(Options,'a0')
        mu0 = Options.a0;
        if ~is_finite_real_two_dimensional_numeric_matrix(mu0) || ...
                size(mu0,1)~=nmu || size(mu0,2)~=T
            error(['a0 field of optional options input corresponds to ',...
                'an initial guess at the quadratic programming ',...
                'solution the shadow shocks and so must be a finite ',...
                'real two dimensional numeric matrix with number of ',...
                'rows equal to the number of constraints (different ',...
                'types of shadow shock) and number of cols equal to ',...
                'the number of periods over which constraints might bind'])
        end
    end
    if isfield(Options,'fvalTol')
        fvalTol = Options.fvalTol;
        if ~is_positive_real_numeric_scalar(fvalTol)
            error(['fvalTol field of optional options input must be a ',...
                'scalar for evaluation of whether or not the quadprog ',...
                'solution to the LCP is accceptably close to zero'])
        end
    end
    if isfield(Options,'Quadprog')
        QuadProgOpts = overlay_default_structure(...
            Options.Quadprog,QuadProgOpts);
    end
end

%% SETUP INFO FOR QUADPROG
% The quadratic programming problem setup below returns the vector of
% shocks that minimises 0.5*x'*H*x + f'*x s.t. A*x<=b & x>=0.
H = (M+M');
f = v;
A = -M;
b = v;
lb = zeros(nmu*T,1);

%% VECTORISE INITIAL GUESS
x0 = reshape(mu0',[T*nmu,1]);

%% RUN QUADPROG
QuadProgOpts.ConvexCheck = 'off';
[x,fval,exitFlag] = quadprog(H,f,A,b,[],[],lb,[],x0,QuadProgOpts);

if exitFlag < 1
    error(['The exitflag from quadprog is ' num2str(exitFlag), ...
        '. A non-positive value means that the routine failed.']);
end
if abs(fval) > fvalTol
    error(['The solution returned by quadprog is not admissable ',...
        'because the (absolute) value of the objective (' num2str(fval) ...
        ') is larger than tolerance (' num2str(fvalTol) ').']);
end

%% RETURN RESULT
% mu^{\ast} corresponds to \alpha^{\ast} in the Holden and Paetz paper.
mustar = reshape(x,[T,nmu])';

end

%% FUNCTION TO DEFINE DEFAULT OPTIONS
function DefaultOpts = get_default_options(nmu,T)
% Returns default options for the algorithm.
%
% INPUTS:   
%   -> nmu: number of constraints
%   -> T: number of time periods over which constraints might bind
%
% OUTPUTS
%   -> DefaultOpts: Structure of default options
%       - a0: initial condition for the search (must be specified as an 
%         nmu*T matrix (rows indexing shadow shock, cols periods)
%       - fvalTol: tolerance for accepting the Linear Complementarity 
%         Solution
%       - QuadProg: structure of options for the QP algorithm
%         (see sub-function below for info on valid content)

%% DEFINE DEFAULT GUESS FOR SHADOW ANTICIPATED SHOCKS AS ZEROS
DefaultOpts.mu0 = zeros(nmu,T);

%% DEFINE DEFAULT FOR EVALUATION OF LCP TO BE SATISFIED
% Which requires that the function evaluation on the solutionr returned is
% (approximately) 0.
DefaultOpts.fvalTol = 1e-6;

%% GET DEFAULT OPTIONS FOR QUADPROG
DefaultOpts.QuadProg = get_quadprog_default_options();

end

%% FUNCTION TO DEFINE DEFAULT QUADPROG OPTIONS
function DefaultQuadProgOpts = get_quadprog_default_options()
% Contains a set of default options to run quadprog.
%
% INPUTS:   
%   -> none
%
% OUTPUTS
%   -> DefaultQuadProgOpts: structure of default options

%% DEFINE OPTIONS USING OPTIMSET
% Note that `active-set' option was removed with Matlab R2016b.

DefaultQuadProgOpts = optimset(...
    'MaxIter',10^4,...
    'Display','off');


end