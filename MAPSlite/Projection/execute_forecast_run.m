function NewRunData = execute_forecast_run(Model,RunData,ProvRunData)
% This "super" forecast macro manages the execution of a forecast run.
% It is a single macro that can be used to compute a new forecast run
% regardless of whether the model used is linear state space or non-linear
% backward-looking.
%
% INPUTS:
%   -> Model: valid MAPS model
%   -> RunData: structure of existing forecast run data
%   -> ProvRunData: structure of provisional run data
%
% OUTPUTS:
%   -> NewRunData: structure of new forecast run data
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> unpack_model
%   -> generate_MAPS_exception_add_cause_and_throw
%   -> execute_LSS_model_forecast_run
%   -> execute_NLBL_model_forecast_run
%
% DETAILS:
%   -> This macro is a macro for the forecast macros! It examines the
%      content of the model to determine its type.
%   -> It then calls the appropriate forecast macro to execute the forecast
%      run based on the model type.
%
% NOTES:
%   -> See <> for more details of forecast run execution using linear state
%      space and non-linear backward-looking models.
%
% This version: 21/02/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number & shape of inputs is as expected. All inputs are 
% compulsory.
if nargin < 3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId); 
elseif ~isstruct(RunData)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);     
elseif ~isstruct(ProvRunData)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);  
end

%% DETERMINE MODEL TYPE
% Unpack the model is linear state space field, which is set true for
% linear state space (LSS) models and false for non-linear backward-looking
% (NLBL) models.
try
    modelIsLinearStateSpace = unpack_model(...
        Model,{'modelIsLinearStateSpace'});
catch ModelUnpackE
    errId = ['MAPS:',mfilename,':ClassDeterminationFailure'];
    generate_MAPS_exception_add_cause_and_throw(ModelUnpackE,errId);
end

%% CALL APPROPRIATE FORECAST EXECUTION MACRO
% Call the appropriate forecast execution macro depending on the model
% type.
if modelIsLinearStateSpace
    NewRunData = execute_LSS_model_forecast_run(Model,RunData,ProvRunData);
else
    NewRunData = execute_NLBL_model_forecast_run(...
        Model,RunData,ProvRunData);
end

end