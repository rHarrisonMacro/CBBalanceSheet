function ss = solve_LSS_model_steady_state(SSfunHandle,theta)
% This module evaluates the steady states to produce a numeric vector.
% It evaluates the recurstively reduced system of steady states and
% parameter transformations to produce a numeric vector for the steady
% states.
%
% INPUTS:   
%   -> SSfunHandle: function handle for the recursively reduced steady 
%      state expressions
%   -> theta: numeric vector of parameters
%
% OUTPUTS:  
%   -> ss: numeric vector of steady states & transformed parameters
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> is_finite_real_numeric_column_vector
%   -> generate_MAPS_exception_add_cause_and_throw
%
% DETAILS:  
%   -> This module evaluates a recursively reduced set of steady state &
%      parameter transformations using a numeric vector of model 
%      parameters.
%   -> Its logic relies on the system of equations being recursive and
%      having been reduced such that the system has a vector analytical
%      solution (i.e. the system has been reduced so that it does not need
%      to be evaluated in any particular order).
%
% NOTES:
%   -> See <> for a description of MAPS LSS models, their rules and a
%      description of how they are solved
%   -> This module may be extended in future to allow for non-analytical
%      systems of equations.
%
% This version: 13/05/2011
% Author(s): Alex Haberis and Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of inputs is as expected.
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~strcmp(class(SSfunHandle),'function_handle')
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~is_finite_real_numeric_column_vector(theta)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
end

%% EVALUATE STEADY STATE & PARAMETER TRANSFORMATIONS
% Evaluate the steady state & parameter transformation expressions in the
% function handle using the vector of parameters. Capture any error in
% evaluating these expessions and throw an exception if necessary.
try
    ss = SSfunHandle(theta);
catch ssEvalE
    errId = ['MAPS:',mfilename,':SteadyStateEvalFailure'];
    generate_MAPS_exception_add_cause_and_throw(...
        ssEvalE,errId);
end

end