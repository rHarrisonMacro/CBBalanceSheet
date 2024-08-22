function x = ...
    evaluate_endogenous_vars_in_ES_model_using_linear_interpolation(...
    snodes,f,s)
% Evaluates expectations of the endogenous variables in ES model.
%
% INPUTS:   
%   -> snodes: 1*ds cell array of states across each dimension
%   -> f: ns*dx policy functions matrix of policy function
%   -> s: 1*ds vector of states
%
% OUTPUTS:  
%   -> x: 1*dx vector of endogenous variables

%% GET INDICES AND WEIGHTS TO APPLY
[inds,weights] = lookup_indices_and_weights_for_linear_interpolation(...
    s,snodes);

%% COMPUTE ENDOGENOUS VARIABLES VIA LINEAR INTERPOLATION & EXTRAPOLATION
x = weights'*f(inds,:);

end