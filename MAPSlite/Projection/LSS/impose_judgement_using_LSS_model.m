function NewRunData = impose_judgement_using_LSS_model(...
    Model,RunData,Judgements,algorithmIsMinimumDamage)
% This macro update a LSS model projection for forecast judgements.
% It allows users to impose judgement on a projection using a linear state 
% space (LSS) model via: a) the shocks and/or b) by fixing endogenous 
% variables (in which case each fix can be treated as either anticipated 
% or unanticipated) and/or c) by changing any time-varying trends in the 
% model.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> RunData: existing forecast run dataset structure
%   -> Judgements: structure with info about all judgements to be imposed
%       - AnticipatedFixes (optional): structure of anticipated fixes
%           - modelVariables: a cell array of metadata & fix values
%           - modelObservables: a cell array of metadata & fix values
%           - rawObservables: a cell array of metadata & fix values
%           - shockUsages: a cell array of metadata & usage indicators
%       - UnanticipatedFixes (optional): structure of unanticipated fixes
%           - modelVariables: a cell array of metadata & fix values
%           - modelObservables: a cell array of metadata & fix values
%           - rawObservables: a cell array of metadata & fix values
%           - shockUsages: a cell array of metadata & usage indicators
%       - Shocks (optional): structure detailing judgements over shocks:
%           - anticipated: a cell array of metadata & judgements
%           - unanticipated: a cell array of metadata & judgements
%       - timeVaryingTrends (optional): cell array of metadata & judgements
%   -> algorithmIsMinimumDamage (optional): true or false for the inversion
%      algorithm used
%
% OUTPUTS:  
%   -> NewRunData: structure describing new judgemental forecast
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
%   -> validate_run_dataset_for_LSS_model_project_macros
%   -> validate_LSS_model_judgements_dataset
%   -> update_time_trends_forecast_dataset_with_judgements
%   -> update_shocks_forecast_dataset_with_judgements
%   -> unpack_model
%   -> create_LSS_model_inversion_instructions
%   -> invert_LSS_model
%   -> project_using_LSS_model
%   -> validate_forecast_run_dataset_against_judgements (to be added)
%
% DETAILS:  
%   -> This linear state space model forecast macro manages the imposition 
%      of five diffierent types of forecast judgement: (i) judgements to 
%      raw observable; (ii) judgements to model observables; (iii) 
%      judgements to model variables; (iv) judgements to the underlying 
%      shocks; (v) judgements to deterministic time-varying trends.
%   -> Judgements to the shocks (iv) or any time-varying trends (v)
%      are applied directly to the matrices of existing, pre-judgemental 
%      shocks and trends. 
%   -> Judgements made to any of the endogenous variables (i)-(iii) can be 
%      applied either as anticipated (if the model is forward looking) or
%      unanticipated. All of these judgements are combined in model
%      variables space and then the model is inverted to find a new set of
%      of forecast shocks which suuport the projection.
%   -> This new set of forecast shocks (which incorporates any judgements 
%      made in shock space) is combined with the input constraints
%      (forecast initial conditions) to update the projection (with any
%      time-varying trend judgements made affecting the projection for raw
%      observable variables).
%   -> This macro throws an exception if the inputs are not as expected, if 
%      they are not consistent with each other or if any of the LSS model 
%      impose judgement rules are broken. 
%
% NOTES:
%   -> See the MAPS technical documentation for more details on this macro
%      and the modules called herewith.
%   -> The algebra and rules for impose judgement are documented in
%      ANALYTCAL xxxxxxx.
%
% This version: 04/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS AGAINST EXPECTATIONS
% Complete basic checks on the number and shapes of inputs. Throw an error
% if they are not as expected.
if nargin < 3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(RunData)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(Judgements)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
elseif nargin>3 && (~isscalar(algorithmIsMinimumDamage)||...
        ~islogical(algorithmIsMinimumDamage))
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);    
end

%% VALIDATE MODEL & DATASET PASSED IN
% Validate dataset passed in for use in this projection macro. It will
% throw an error if the content of the dataset is not consistent with the 
% requirements of this macro or the model being used.
validate_run_dataset_for_LSS_model_project_macros(Model,RunData);

%% VALIDATE JUDGEMENT DATASET PASSED IN
% Valdiate the judgements dataset structure passed in to check that the 
% shape and content of the judgements is as expected and is consistent with 
% the model input. And that the judgements being made satisfy the rules for 
% imposing judgements. If they do not, this function throws an exception 
% detailing which rule(s) were broken and why.
[JudgementsInclude,JudgementsInfo] = ...
    validate_LSS_model_judgements_dataset(Model,Judgements);

%% SETUP OUTPUT
% Set the output to be updated below as equal to the data input.
NewRunData = RunData;

%% UPDATE TIME-VARYING TRENDS WITH JUDGEMENTS (IF APPLICABLE)
% Update the existing, pre-judgemental deterministic time-varying trend 
% forecast data with any time-varying trend judgements being made.
if JudgementsInclude.timeVaryingTrendJudgements
    NewRunData.Forecast.timeVaryingTrends = ...
        update_time_trends_forecast_dataset_with_judgements(Model,...
        RunData.Forecast.timeVaryingTrends,Judgements.timeVaryingTrends);
end

%% UPDATE SHOCKS WITH JUDGEMENTS (IF APPLICABLE)
% Update the existing, pre-judgemental shocks forecast data with any shock 
% judgements being made.
if JudgementsInclude.anticipatedShockJudgements || ...
        JudgementsInclude.unanticipatedShockJudgements
    NewRunData.Forecast.Shocks = ...
        update_shocks_forecast_dataset_with_judgements(...
        Model,RunData.Forecast.Shocks,Judgements.Shocks);
end

%% UPDATE SHOCKS USED TO IMPLEMENT FIXES (IF APPLICABLE)
% If any fix judgements were input, invert the linear state space model to 
% find new forecast values for the shocks used in the inversion. Unpack the
% components of the model and the forecast run required for the LSS model
% inversion. Create the inversion instructions from judgements dataset
% structure and, if not provided as optional input to this function, set 
% the indicator for the algorithm to use in the inversion to be "minimum 
% damage" (see MAPS documentation for details). Call the LSS model 
% inversion module to update the forecast values for the shocks and then 
% pack them into the new forecast run dataset.
if JudgementsInclude.shocksForAnticipatedFixes || ...
        JudgementsInclude.shocksForUnanticipatedFixes
    [B,PHI,F] = unpack_model(Model,{'B','PHI','F'});
    xT = NewRunData.Constraint.modelVariables;
    Shocks = NewRunData.Forecast.Shocks;
    InversionInstructions = create_LSS_model_inversion_instructions(...
        Model,NewRunData,JudgementsInclude,JudgementsInfo);
    if nargin < 4
        InversionInstructions.algorithmIsMinimumDamage = true;
    end
    NewShocks = invert_LSS_model(B,PHI,F,xT,Shocks,InversionInstructions);
    NewRunData.Forecast.Shocks = NewShocks;    
end

%% CALL THE PROJECTION MACRO TO UPDATE THE PROJECTION FOR NEW SHOCKS
% Update the projections in the new run for the judgements made (as
% reflected in new shock values and/or new time-varying trends) using the
% MAPS LSS model project macro.
NewRunData = project_using_LSS_model(Model,NewRunData);

%% COMPARE THE NEW RUN DATA TO FIXES (IF APPLICABLE)
% If any fix judgements were input, compare their values with the outcome
% from the inversion and projection. Issue warnings if any of the fixes
% failed to hold within tolerance.
if isfield(Judgements,'AnticipatedFixes') || ...
        isfield(Judgements,'UnanticipatedFixes')
%     validate_forecast_run_dataset_against_judgements(...
%         NewRunData,Judgements);
end

end