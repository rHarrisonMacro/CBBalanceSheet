function [ShocksInclude,H,nz] = interrogate_forecast_shocks_structure(...
    Shocks,PHI,F)
% This helper interrogates a LSS model forecast shocks structure. 
% It examines the content of a structure containing anticipated and/or
% unanticipated shocks to valdiate whether or not they are valid for use in
% MAPS linear state space (LSS) model forecast modules. It can also be used
% to return information about the content of the shocks and the forecast
% horizon that can be inferred from them.
%
% INPUTS:
%   -> Shocks: struture with: 
%       - anticipated (model dependent): nz*H matrix of anticipated shocks
%       - unanticipated: nz*H matrix of unanticipated shocks
%   -> PHI: matrix of loadings on shocks from the model solution
%   -> F: matrix of loadings on anticipated shocks from the model solution
%
% OUTPUTS
%   -> ShocksInclude: struture with: 
%       - anticipated: true if input contains anticipated shocks, false
%         otherwise
%       - unanticipated: true if input contains unanticipated shocks, false
%         otherwise
%   -> H: numeric scalar representing the forecast horizon 
%   -> nz: dimensions of the shocks
%
% CALLS: 
%   -> generate_and_throw_MAPS_exception
%   -> is_finite_real_two_dimensional_numeric_matrix
%
% DETAILS: 
%   -> This helper validates a set of shocks for use in MAPS forecast
%      modules. It checks the content of the structure input and checks
%      consistency of the shocks with each other.
%
% NOTES:   
%   -> See <> for a description of MAPS forecast helpers and data 
%      validation.
%
% This version: 14/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of the inputs is as expected.
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});  
elseif ~isstruct(Shocks)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
end

%% COMPUTE DIMENSIONS OF THE SHOCKS IN THE MODEL
% Compute the dimensions of the shocks in the model as the number of 
% columns in the shock loadings matrix from the solution, where:
% x{t} = B*x{t-1}+PHI*z{t}  
nz = size(PHI,2);

%% CHECK SHOCKS
% Check for the existence of each shock type in the input. If it does exist
% check that the data is in a finite, numeric two-dimensional matrix.
shockTypes = {'anticipated';'unanticipated'};
nShockTypes = size(shockTypes,1);
for iType = 1:nShockTypes
    iShockType = shockTypes{iType};
    if isfield(Shocks,iShockType)
        ShocksInclude.(iShockType) = true;
        if ~is_finite_real_two_dimensional_numeric_matrix(Shocks.(iShockType))
            errId = ['MAPS:',mfilename,':BadShocksData'];
            generate_and_throw_MAPS_exception(errId,{iShockType});
        elseif size(Shocks.(iShockType),1) ~= nz
            errId = ['MAPS:',mfilename,':IncompatiblePHIShocks'];
            generate_and_throw_MAPS_exception(errId,{iShockType});
        end
    else
        ShocksInclude.(iShockType) = false;
    end
end

%% CHECK CONSISTENCY OF SHOCKS WITH EACH OTHER
% If the shocks input contains both anticipated and unanticipated shocks
% then check that the dimensions are consistent with each other. Throw an
% exception if not, or if anticipated shocks are included but the model is
% backward-looking, or if neither type of shock exists.
if ShocksInclude.anticipated
    H = size(Shocks.anticipated,2);
    if ShocksInclude.unanticipated
        if size(Shocks.unanticipated,2) ~= H
            errId = ['MAPS:',mfilename,':InconsistentShocksTimeDims'];
            generate_and_throw_MAPS_exception(errId);
        end
    elseif ~any(any(F))
        errId = ['MAPS:',mfilename,':BadUseOfAnticipatedShocks'];
        generate_and_throw_MAPS_exception(errId);        
    end
elseif ShocksInclude.unanticipated
    H = size(Shocks.unanticipated,2);
else
    errId = ['MAPS:',mfilename,':MissingShocks'];
    generate_and_throw_MAPS_exception(errId);
end

end