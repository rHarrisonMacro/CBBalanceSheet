function [minStateInds,minStateTransformMat] = ...
    compute_minimum_state_space_representation(Model)
% This function computes the minimum state representation of an LSS model.
% It returns index numbers and a transformation matrix that can be used to
% transform the input model into its minimum state space equivalent. This
% can be used to improve the performance of various model operations. 
%
% INPUTS:   
%   -> Model: MAPS LSS model object
%
% OUTPUTS:  
%   -> minStateInds: indices of the model variables in min state
%   -> minStateTransformMat: matrix equivalent of the indices
%
% DETAILS:  
%   -> This function returns all the information necessary to compute the 
%      minimum state space representation of an LSS model.
%   -> It works by dropping variables that do not appear in model equations
%      with a lag and that do not feature in the measurement equations (if
%      applicable).
%   -> The set of minimum state model mnemonics can then be computed as
%      xMinStateMnems = xMnems(minStateInds).
%   -> And model solution or structural matrices can be transformed as
%      BminState = minStateTransformMat*B*minStateTransformMat'
%           
% NOTES:    
%   -> This function is used in areas of MAPS where performance is
%      important. For example, it is used in the MAPS estimation toolkit to
%      speed-up Kalman filter evaluations in posterior optimisation and
%      simulation.
%           
% This version: 12/12/2013
% Author(s): Richard Harrison, Francesca Monti & Matt Waldron

%% CHECK INPUTS
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
end

%% CHECK IF MODEL HAS MEASUREMENT EQUATIONS
modelHasMeasurementEqs = unpack_model(Model,'modelHasMeasurementEqs');

%% CREATE SYMBOLIC MATRICES
% This ensures that the checks below are analytical rather than numeric.
if modelHasMeasurementEqs
    [HBsym,~,~,~,~,Gsym] = create_LSS_structural_symbolic_matrices(Model);
else
    HBsym = create_LSS_structural_symbolic_matrices(Model);
end

%% FIND MINIMUM STATE SPACE INDICES
% Variables form part of the minimum state space representation if they
% appear with a lag in the model equations or if they appear in the
% measurement equations (if applicable).
minStateLogicals = any(HBsym~=0)';
if modelHasMeasurementEqs
    minStateLogicals = (minStateLogicals|any(Gsym~=0)');
end
minStateInds = find(minStateLogicals);

%% COMPUTE MINIMUM STATE SPACE TRANSFORMATION MATRIX
% The transformation matrix equivalent to the indices is the rows of the
% transformation matrix that correspond to those indices.
nx = size(HBsym,2); 
nxeye = eye(nx);
minStateTransformMat = nxeye(minStateLogicals,:);

end