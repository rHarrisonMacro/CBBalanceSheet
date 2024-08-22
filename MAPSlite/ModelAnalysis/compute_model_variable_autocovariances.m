function [xxCellArray,xxMat] = ...
    compute_model_variable_autocovariances(Model,H,computeAsCorrelation)
% This model analysis function computes covariances of model variables.
% It computes the autocovariances of the variables in the model up to a 
% chosen horizon.
% 
% INPUTS:
%   -> Model: LSS model structure
%   -> H: horizon up to which to compute autocovariances
%   -> computeAsCorrelation: false/true
%
% OUTPUTS:
%   -> xxCellArray: (nx*nx)*3 cell array with model variable mnemonic
%      string pairs in the first two columns and a vectors of 
%      autocovariances between them from the model in the third column
%   -> xxMat: nx*nx*(H+1) matrix of variable autocovarainces from model
%
% DETAILS:
%   -> This model analysis function computes the autocovariances of a 
%      model's variables.
%   -> Using the standard MAPS notation, the autocovariances of the model
%      variables can be defined as E(xi{t},xj{t-h})
%   -> The outputs are such that the first column of the cell array
%      contains the time t dated variable from the definition above (and
%      this is referenced by index in the output matrices), while the 
%      second column contains the lagged variables from the definition 
%      above (with analogous indexing into the matrices).
%   -> The third input is optional and can be used to specifiy the outputs
%      as correlations rather than covariances. The default is for
%      covariances to be returned, so computeAsCorrelation is false.
%   -> This function first checks the inputs. It then computes the 
%      covariances from the model, before packing the cell array for 
%      output.
%
% NOTES:
%   -> Please see the MAPS user guide for more information on model
%      analysis functionality in MAPS.
%
% This version: 24/06/2013
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~is_non_negative_real_integer(H)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif nargin>2 && ~is_logical_scalar(computeAsCorrelation)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);        
end

%% HANDLE OPTIONAL INPUT
if nargin < 3
    computeAsCorrelation = false;
end

%% CHECK MODEL IS COMPATIBLE WITH CALL TO THIS FUNCTION
modelIsLinearStateSpace = unpack_model(Model,{'modelIsLinearStateSpace'});
if ~modelIsLinearStateSpace
    errId = ['MAPS:',mfilename,':BadModelClass'];
    generate_and_throw_MAPS_exception(errId);
end

%% UNPACK REQUIRED MODEL COMPONENTS
[B,P,xMnems] = unpack_model(Model,{'B','P','xMnems'});

%% INITIALISE OUTPUT
nx = size(B,1);
xxMat = NaN*ones(nx,nx,H+1);
xxCellArray = cell(nx*nx,3);

%% COMPUTE AUTO-COVARIANCES IN MODEL
for h = 1:H+1
    xxMat(:,:,h) = B^(h-1)*P;
    if computeAsCorrelation
        if h==1
            xStdVec = diag(P).^0.5;
            xxStdScalingMat = xStdVec*xStdVec';
        end
        xxMat(:,:,h) = xxMat(:,:,h)./xxStdScalingMat;
    end
end

%% PACK AUTOCOVARIANCES INTO CELL
for ix1 = 1:nx
    for ix2 = 1:nx
        xxCellArray((ix1-1)*nx+ix2,:) = ...
            [xMnems(ix2) xMnems(ix1) ...
            {reshape(xxMat(ix1,ix2,:),[1 H+1])}];
    end
end

end