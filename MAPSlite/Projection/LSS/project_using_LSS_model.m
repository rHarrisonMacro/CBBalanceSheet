function NewRunData = project_using_LSS_model(Model,RunData)
% This macro produces a projection for all variables in a LSS model.
% This forecast macro computes or updates a projection using a linear state 
% space (LSS) model given a set of forecast initial conditions 
% (constraints) and a set of existing forecast judgements (as summarised by 
% the values in the model shocks over the forecast).
%
% INPUTS:
%   -> Model: MAPS LSS model structure
%   -> RunData: structure of forecast run data with the following fields
%       - Constraint: structure of forecast initial conditions
%           - modelVariables: vector of model variable constraints
%           - rawObservables (model dep.): vector of raw observable 
%             constraints
%           - timeVaryingTrends (model dep.): vector of time-varying trend 
%             constraints
%       - Forecast: stucture of forecast data with the following fields
%           - Shocks: structure of existing shocks
%               - unanticipated: matrix of unanticipated forecast shocks
%               - anticipated (model dep.): matrix of anticipated
%                 forecast shocks
%           - timeVaryingTrends (model dep.): matrix of time-varying trends
%             over the forecast
%
% OUTPUTS:
%   -> NewRunData: structure of new forecast run data
%       - Constraint: structure of forecast initial conditions
%           - modelVariables: vector of model variable constraints
%           - rawObservables (model dep.): vector of raw observable 
%             constraints
%           - timeVaryingTrends (model dep.): vector of time-varying trend 
%             constraints
%       - Forecast: stucture of forecast data with the following fields
%           - modelVariables: matrix of forecasts for the model variables
%           - modelObservables (model dep.): matrix of forecasts for the 
%             model observables
%           - rawObservables: matrix of forecasts for the raw observables
%           - Shocks: structure of existing shocks
%               - unanticipated: matrix of unanticipated forecast shocks
%               - anticipated (model dep.): matrix of anticipated
%                 forecast shocks
%           - timeVaryingTrends (model dep.): matrix of time-varying trends
%             over the forecast
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> unpack_model
%   -> validate_run_dataset_for_LSS_model_project_macros
%   -> project_LSS_model_variables
%   -> compute_observables_from_model_variables
%   -> transform_observables_from_model_to_raw_space
%   -> validate_retransformed_raw_observable_dataset
%
% DETAILS:
%   -> This macro produces a projection for all variables in the linear 
%      state space model input.
%   -> It first computes a projection for the model variables given initial
%      conditions and a set of shocks over the forecast.
%   -> If the model has measurement equations, it then computes projections
%      for the model observables consistent with the model variables.
%   -> Finally and if the model has data transformation equations, it
%      computes projections for the raw observables consistent with the 
%      model observables.
%
% NOTES:
%   -> See <> for more details of forecast run execution using linear state
%      space models.
%   -> Note that if the run data represents a complete, consistent existing
%      forecast run then this macro will just replicate that forecast run.
%
% This version: 02/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number & shape of inputs is as expected. All inputs are 
% compulsory.
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId); 
elseif ~isstruct(RunData)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);   
end
    
%% VALIDATE DATASET PASSED IN
% Validate dataset passed in for use in this projection macro. It will
% throw an error if the content of the dataset is not consistent with the 
% requirements of this macro or the model being used. 
validate_run_dataset_for_LSS_model_project_macros(Model,RunData);

%% UNPACK MODEL CHARACTERISTICS INFO
% Unpack model characteristics information (i.e. whether or not the LSS 
% model has measurement & data transformation equations).
[modelHasMeasurementEqs,modelHasDataTransformationEqs] = unpack_model(...
    Model,{'modelHasMeasurementEqs','modelHasDataTransformationEqs'});

%% SETUP OUTPUT
% Set the output to be updated below as equal to the data input.
NewRunData = RunData;

%% PROJECT MODEL VARIABLES
% Unpack the model solution matrices, the model variable initial conditions 
% (projection constraints) and the forecast shocks. Call the model variable
% projection module and pack the result into the output.
[B,PHI,F] = unpack_model(Model,{'B','PHI','F'});
xT = RunData.Constraint.modelVariables;
Shocks = RunData.Forecast.Shocks;
xf = project_LSS_model_variables(B,PHI,F,xT,Shocks);
NewRunData.Forecast.modelVariables = xf;

%% PROJECT MODEL OBSERVABLES
% If the model has measurement equations unpack the measurement equation 
% matrices and compute the model observables consistent with the model 
% variable projection from above.
if modelHasMeasurementEqs
    [D,G] = unpack_model(Model,{'D','G'});
    Yf = compute_observables_from_model_variables(xf,D,G);
    NewRunData.Forecast.modelObservables = Yf;
end

%% PROJECT RAW OBSERVABLES
% If the model has data transformation equations unpack the data
% retransformation function handle and the time trend characteristics
% field. Unpack the raw observable initial conditions and then retransform
% the observables from model to raw space using the time trends in the data 
% input if they are part of the model input and without them if not. Pack
% the resulting raw observable forecasts into the data output.
if modelHasDataTransformationEqs
    [RTfunHandle,modelHasTimeVaryingTrends] = unpack_model(...
        Model,{'RTfunHandle','modelHasTimeVaryingTrends'});
    YtildeT = RunData.Constraint.rawObservables;
    if modelHasTimeVaryingTrends
        etatT = RunData.Constraint.timeVaryingTrends;
        etatf = RunData.Forecast.timeVaryingTrends;
        Ytildef = transform_observables_from_model_to_raw_space(...
            RTfunHandle,Yf,YtildeT,etatf,etatT);
    else
        Ytildef = transform_observables_from_model_to_raw_space(...
            RTfunHandle,Yf,YtildeT);  
    end
    validate_retransformed_raw_observable_dataset(Model,Ytildef);
    NewRunData.Forecast.rawObservables = Ytildef;
end

end