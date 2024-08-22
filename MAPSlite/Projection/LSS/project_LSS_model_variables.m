function xf = project_LSS_model_variables(B,PHI,F,xT,Shocks)
% This module produces a projection for LSS model variables.
% The projection is a judgemental projection (where the judgements are 
% summarised by a set of unaticipated and possibly anticipated shocks) 
% conditional on a set of initial conditions using a MAPS linear state 
% space (LSS) model.
%
% INPUTS:
%   -> B: nx*nx matrix of loadings on lagged model variables
%   -> PHI: nx*nz matrix of loadings on shocks
%   -> F: nx*nx matrix of loadings on anticipated shocks
%   -> xT: nx*1 vector of initial conditions for model variables
%   -> Shocks: struture with: 
%       - anticipated (model dependent/optional): nz*H matrix of shocks
%       - unanticipated (optional): nz*H matrix of shocks
%
% OUTPUTS
%   -> xf: nx*H matrix of projected model variables
%
% CALLS: 
%   -> generate_and_throw_MAPS_exception
%   -> validate_LSS_model_solution_matrices
%   -> interrogate_LSS_model_variable_initial_conditions
%   -> interrogate_forecast_shocks_structure
%   -> compute_anticipated_shocks_impact (sub-function)
%
% DETAILS: 
%   -> This module projects LSS model variables forward given a set of 
%      initial conditions and a set of judgements summarised in values for
%      unanticipated and/or anticipated shocks.
%   -> It caters for three cases. First, the input contains forecast data 
%      for unanticipated shocks only (perhaps because the model is 
%      backward-looking). Second, the input contains forecast data for 
%      anticipated shock only. Third, the input contains forecast data for 
%      both types of shock (the most common case for judgemental 
%      projections with forward-looking LSS models).
%   -> The forecast horizon is inferred from the number of data points in
%      the shock matrices (which must be consistent for anticipated and
%      unanticipated shocks if both are passed in).
%   -> For a projection method without any judgement (i.e. no shocks), see 
%      plain_vanilla_project_LSS_model_variables.
%
% NOTES:   
%   -> See <> for a description of MAPS forecast modules and macros.
%
% This version: 24/02/2011
% Author(s): Francesca Monti & Matt Waldron

%% CHECK INPUTS
% Check that the number of inputs is as expected. All inputs are 
% compulsory.
if nargin < 5
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
end

%% VALIDATE SOLUTION MATRICES
% Call a MAPS module helper to validate that the solution matrices are as 
% expected and consistent with one another.
validate_LSS_model_solution_matrices(B,PHI,F);

%% VALIDATE MODEL VARIABLE INITIAL CONDITIONS
% Call a MAPS module helper to validate that the projection initial 
% conditions are as expected and consistent with the model solution. This
% helper returns the dimension of the model variables.
nx = interrogate_LSS_model_variable_initial_conditions(xT,B);

%% VALIDATE SHOCKS INPUT
% Call a forecast module helper to validate the shocks input. This function
% returns a structure describing the content of the shocks structure and a
% numeric scalar representing the forecast horizon. It will throw an error 
% if the dimensions of the shocks are not consistent with each other, if
% they are not consistent with the model solution or if anticipated shocks
% were input but the model is backward looking ("F" is all zeros).
[ShocksInclude,H] = interrogate_forecast_shocks_structure(Shocks,PHI,F);
        
%% PROJECT MODEL VARIABLES
% Project the model variables using the model plus judgement on anticipated 
% and/or unanticipated shocks depending on what was input.  Note also that  
% the projection maths is repeated across each case of the three valid 
% cases. In this instance, this repetition is desirbale because it imporves 
% performance.
xf = zeros(nx,H);
if ~ShocksInclude.anticipated && ShocksInclude.unanticipated
    u = Shocks.unanticipated;
    uImpact = zeros(nx,1);
    for s = 1:H
        uImpact = B*uImpact+PHI*u(:,s);
        xf(:,s) = B^s*xT+uImpact;
    end
elseif ShocksInclude.anticipated && ~ShocksInclude.unanticipated
    a = Shocks.anticipated;
    R = compute_anticipated_shocks_impact(PHI,F,a,H,nx);
    aImpact = zeros(nx,1);
    for s = 1:H
        aImpact = B*aImpact+R(:,s);
        xf(:,s) = B^s*xT+aImpact;
    end
else
    a = Shocks.anticipated;
    u = Shocks.unanticipated;
    R = compute_anticipated_shocks_impact(PHI,F,a,H,nx);
    aImpact = zeros(nx,1);
    uImpact = zeros(nx,1);
    for s = 1:H
        aImpact = B*aImpact+R(:,s);
        uImpact = B*uImpact+PHI*u(:,s);
        xf(:,s) = B^s*xT+aImpact+uImpact;
    end  
end
    
end

%% FUNCTION TO COMPUTE ANTICIPATED SHOCKS IMPACT
function R = compute_anticipated_shocks_impact(PHI,F,a,H,nx)
% This helper computes the impact of the anticipated shocks. 
% Specifically, it computes the R vector for each period of the projection.
% For each period in the projection (given by each column of the R matrix),
% the R vector describes the impact of all of the anticipated shocks in the
% projection on the model variables projection.
%
% INPUTS:   
%   -> PHI: loadings on the shocks in the model solution
%   -> F: forward loadings to "discount" anticipated shocks
%   -> a: matrix of anticipated shocks over the projection
%   -> H: horizon for the projection
%   -> nx: dimensions of the model variables
%
% OUTPUTS:  
%   -> R: marix of impacts of all the anticipated shocks in each period
%
% CALLS:
%   -> none

%% INITIALISE THE OUTPUT AND THE LATEST VALUES FOR THE IMPACTS
R0 = zeros(nx,1);
R = zeros(nx,H);

%% COMPUTE R BY BACKWARD ITERATION
for s = H:-1:1
    R(:,s) = F*R0+PHI*a(:,s);
    R0 = R(:,s);    
end

end