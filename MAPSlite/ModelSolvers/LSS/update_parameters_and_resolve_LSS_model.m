function UpdatedModel = update_parameters_and_resolve_LSS_model(...
    Model,thetaUpdate,thetaUpdateIdentifiers,theta,thetaMnems)
% This LSS model helper updates the parameters of an LSS model & resolves. 
% It can be operated in a number of different modes depending on the 
% information passed in, allowing this function to operate at higher and 
% lower levels. 
% 
% INPUTS:
%   -> Model: MAPS LSS model structure
%   -> thetaUpdate: values for the parameters to be updated
%   -> thetaUpdateIdentifiers: identifiers for the parameters to be updated
%      (either index numbers of mnemonics)
%   -> theta (optional): column vector of all parameter values in model
%   -> thetaMnems (optional): column cell string array of all paarameter 
%      mnemonics in model
%
% OUTPUTS:
%   -> UpdatedModel: MAPS LSS model structure resolved with updated params
%
% DETAILS:
%   -> This LSS model utility function updates a subset of the model's
%      parameters and resolves the model.
%   -> The parameters to update can be identified by either index numbers
%      (which will be quickest) or mnemonics in which case the index
%      numbers are looked up.
%   -> In addition, the full vector of parameters may be passed in (which
%      saves unpacking them in the update) and the full set of parameter
%      mnemonics which will also save unpacking them to compute the index
%      numbers (if necessary).
%
% NOTES:
%   -> This helper is useful in MAPS LSS model estimation.
%
% This version: 03/12/2012
% Author(s): Matt Waldron

%% CHECK INPUTS
% The input checking is left to the function calls below.
if nargin < 3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
end

%% UPDATE PARAMETER VECTOR
if nargin < 4
    thetaUpdated = update_parameters(...
        Model,thetaUpdate,thetaUpdateIdentifiers);
elseif nargin < 5
    thetaUpdated = update_parameters(...
        Model,thetaUpdate,thetaUpdateIdentifiers,theta);
else
    thetaUpdated = update_parameters(...
        Model,thetaUpdate,thetaUpdateIdentifiers,theta,thetaMnems);
end

%% RESOLVE MODEL
try
    UpdatedModel = solve_LSS_model(Model,thetaUpdated);
catch SolveE
    errId = ['MAPS:',mfilename,':UnableToResolve'];
    generate_MAPS_exception_add_cause_and_throw(SolveE,errId);
end

end  