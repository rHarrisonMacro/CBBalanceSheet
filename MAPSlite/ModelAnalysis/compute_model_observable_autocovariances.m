function [YYcellArray,YYmodel,YYdata] = ...
    compute_model_observable_autocovariances(Model,Y,H)
% This model analysis function computes covariances of model observables.
% It computes the autocovariances of the observables in the model up to a 
% chosen horizon and computes the equivalent covariances in the data.
% 
% INPUTS:
%   -> Model: LSS model structure
%   -> Y: nY*T matrix of model observables over a particular sample
%   -> H: horizon up to which to compute autocovariances
%
% OUTPUTS:
%   -> YYcellArray: (nY*nY)*4 cell array with model observable mnemonic
%      string pairs in the first two columns and a vectors of 
%      autocovariances between them from the model and data in the third 
%      and fourth columns respectively
%   -> YYmodel: nY*nY*(H+1) matrix of observable autocovarainces from model
%   -> YYdata: nY*nY*(H+1) matrix of observable autocovarainces from data
%
% DETAILS:
%   -> This model analysis function computes the autocovariances of a 
%      model's observables and the equivalent set of covariances implied by
%      the data.
%   -> Using the standard MAPS notation, the autocovariances of the model
%      observables can be defined as E(Yi{t},Yj{t-h})
%   -> The outputs are such that the first column of the cell array
%      contains the time t dated variable from the definition above (and
%      this is referenced by index in the output matrices), while the 
%      second column contains the lagged variables from the definition 
%      above (with analogous indexing into the matrices).
%   -> This function first checks the inputs. It then computes the 
%      covariances from the model and the data in turn, before packing the
%      cell array for output.
%
% NOTES:
%   -> Please see the MAPS user guide for more information on model
%      analysis functionality in MAPS.
%
% This version: 16/03/2013
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~is_finite_real_two_dimensional_numeric_matrix(Y)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif ~is_non_negative_real_integer(H)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);    
end

%% CHECK MODEL IS COMPATIBLE WITH CALL TO THIS FUNCTION
modelIsLinearStateSpace = unpack_model(Model,{'modelIsLinearStateSpace'});
if ~modelIsLinearStateSpace
    errId = ['MAPS:',mfilename,':BadModelClass'];
    generate_and_throw_MAPS_exception(errId);
end
modelHasMeasurementEqs = unpack_model(Model,{'modelHasMeasurementEqs'});
if ~modelHasMeasurementEqs
    errId = ['MAPS:',mfilename,':BadModelCharacteristic'];
    generate_and_throw_MAPS_exception(errId);
end

%% UNPACK REQUIRED MODEL COMPONENTS
% Unpack all components of the model required for what follows including,
% if applicable, the measurement errors.
[B,P,G,Ymnems,modelHasMeasurementErrors] = unpack_model(...
    Model,{'B','P','G','Ymnems','modelHasMeasurementErrors'});
if modelHasMeasurementErrors
    V = unpack_model(Model,'V');
end

%% CHECK DIMENSIONS OF OBSERVABLES & DATASET
[nY,T] = size(Y);
if size(G,1) ~= nY
    errId = ['MAPS:',mfilename,':DataIncompatibleWithModel'];
    generate_and_throw_MAPS_exception(errId);
end
if T < H+1
    errId = ['MAPS:',mfilename,':NotEnoughTimePeriodsInData'];
    generate_and_throw_MAPS_exception(errId);
end

%% INITIALISE OUTPUT
YYmodel = NaN*ones(nY,nY,H+1);
YYdata = NaN*ones(nY,nY,H+1);
YYcellArray = cell(nY*nY,4);

%% COMPUTE AUTO-COVARIANCES IN MODEL
for h = 1:H+1
    YYmodel(:,:,h) = G*B^(h-1)*P*G';
    if h==1 && modelHasMeasurementErrors
        YYmodel(:,:,1) = YYmodel(:,:,1)+V*V';
    end
end

%% COMPUTE AUTO-COVARIANCES IN DATA
for h = 1:H+1
    YYdata(:,:,h) = (1/(T-h))*(Y(:,h:T)-mean(Y(:,h:T),2)*ones(1,T-h+1))*...
        (Y(:,1:T-h+1)-mean(Y(:,1:T-h+1),2)*ones(1,T-h+1))';
end

%% PACK AUTOCOVARIANCES INTO CELL
for iY1 = 1:nY
    for iY2 = 1:nY
        YYcellArray((iY1-1)*nY+iY2,:) = ...
            [Ymnems(iY2) Ymnems(iY1) ...
            {reshape(YYmodel(iY1,iY2,:),[1 H+1])} ...
            {reshape(YYdata(iY1,iY2,:),[1 H+1])}];
    end
end

end