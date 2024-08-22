function [fSlice,sSlice] = slice_policy_function(f,s,sDim,sVal)
%slice_policy_function This code produces a slice of a policy function,
%holding the value of one of the states fixed.
%   
% INPUTS
% -> f: ns*dx matrix of policy functions
% -> s: ns*ds matrix of grid values for the states
% -> sDim: scalar (between 1 and ds inclusive) index of state variable to
% be fixed.
% -> sVal: value of state variable in grid to fix to
%
% OUTPUTS
% -> fslice: nstilde*dx matrix of policy function slice (nstilde is
% ns/number of nodes in sDim)
% -> sSlice: nstilde*(ds-1) matrix of grid values for remaining states

%% INTERROGATE INPUTS
[~,ds] = size(s);
if sDim<0 || sDim>ds
    error(['State dimension to condition on must lie between 1 and ',...
        'total number of states']);
end
sToTest = s(:,sDim);

%% COMPUTE INDICES OF ROWS TO EXTRACT
if sum(sToTest==sVal)==0
    error(['Value to condition on must be a grid point for the ',...
        'selected state']);
end

%% EXTRACT RELEVANT SLICE
sliceIndices = (sToTest==sVal);
fSlice = f(sliceIndices,:);
sSlice = s(sliceIndices,:);
sSlice(:,sDim) = [];

end

