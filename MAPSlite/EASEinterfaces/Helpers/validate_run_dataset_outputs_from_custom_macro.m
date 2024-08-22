function validate_run_dataset_outputs_from_custom_macro(...
    Model,runDatasets,provRunDatasets,BaseRunData,BaseProvRunData)
% This helper validates the outputs of a custom macro.
% It validates both the forecast run dataset and the provisional run
% dataset outputs from a custom macro to ensure that invalid outputs are
% never passed back to EASE.
%
% INPUTS:
%   -> Model: MAPS model structure
%   -> runDatasets: nRuns*1 cell array of forecast run outputs
%   -> provRunDatasets: nRuns*1 cell array of provisional run outputs
%   -> BaseRunData: existing forecast run dataset structure input to macro
%   -> BaseProvRunData: existing provisional run dataset structure input to
%      macro
%
% OUTPUTS:
%   -> none
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> validate_forecast_run_dataset_output
%   -> generate_MAPS_exception_and_add_cause
%   -> validate_provisional_run_dataset_output
%   -> generate_MAPS_exception_and_add_as_cause
%
% DETAILS:
%   -> This helper validates the outputs from a custom macro to protect
%      EASE from invalid data being passed back and to assist MAPS macro
%      builders in building valid macros.
%   -> Custom macros can be used to execute one or more forecast runs.
%      Associated with each forecast run output is a forecast run dataset
%      and a provisional run dataset.
%   -> MAPS forecast run datasets are stored as structures with the fields
%      (and sub fields) of the structures acting as metadata for the
%      series. For example (and if the config information dictates),
%      forecast data for raw observables would be stored as an
%      nYtilde*H matrix (where H is the forecast horizon) in
%      RunData.Forecast.rawObservables.
%   -> MAPS provisional run datasets are stored as structures with the
%      fields (and sub fields) of the structures acting as metadata for the
%      series. For example (and if the config information dictates),
%      an anticipated judgement for a single raw observable would be stored
%      as a 1*2 cell array with the first element containing model metadata
%      for the variable in question (eg 'gdpkp') and the second element
%      containing judgemental data for that variable in a 1*H vector (where
%      H is the forecast horizon) in
%      ProvRunData.Forecast.AnticipatedFixes.rawObservables.
%   -> This function completes a comprehensive validation of both the
%      forecast run and provisional run dataset outputs.
%   -> It first validates that the datasets output are structures in
%      nRuns*1 cell arrays.
%   -> It then calls helper functions to validate the content of each
%      dataset on every run output. The content validation basically checks
%      for two things. First, that the datasets contain the expected
%      variables given the model being used. Second, that the datasets
%      contain valid data (no infs etc) with the correct dimensions.
%   -> Finally, that the datasets are consistent with each other and will
%      not invalidate the treatment of "introduce news" runs in EASE in
%      that if past data is included in a forecast run dataset output, it
%      should also be included in the provisional run dataset output.
%   -> If any of the run datasets output fail validation for any reason,
%      this function will throw an exception detailing the run numbers that
%      failed and the causes of the validation failure.
%
% NOTES:
%   -> See xxxxxx for a discussion of MAPS forecast run datasets and their
%      validation.
%   -> This function leaves all input validation checking to specific
%      checks below and to the functions called.
%
% This version: 21/03/2011
% Author(s): Matt Waldron

%% CHECK NUMBER OF INPUTS
% Check that the number of inpts passed in is as expected by this function.
if nargin < 5
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
end

%% COMPUTE THE TIME PERIOD DIMENSIONS OF THE DATA
% Compute the number of periods in the past and forecast data consistent
% with the base run data input from EASE to MAPS' custom macro interface.
if isfield(BaseProvRunData,'Past')
        T = size(BaseProvRunData.Past.rawObservables,2)-1;
else
        T = size(BaseRunData.Past.modelVariables,2);
end
H = size(BaseRunData.Forecast.Shocks.unanticipated,2);

%% CHECK OUTPUTS
% Check that the outputs from the custom macro are as expected in that they
% are consistent with the assumptions in the cell below which packs the
% EASE outputs.
if ~iscell(runDatasets) || size(runDatasets,2)~=1 || ...
        ndims(runDatasets)~=2 || ~all(cellfun(@isstruct,runDatasets))
    errId = ['MAPS:',mfilename,':BadForecastRunsOutput'];
    generate_and_throw_MAPS_exception(errId);
elseif ~iscell(provRunDatasets) || size(provRunDatasets,2)~=1 || ...
        ndims(provRunDatasets)~=2 || ...
        ~all(cellfun(@isstruct,provRunDatasets))
    errId = ['MAPS:',mfilename,':BadProvRunsOutput'];
    generate_and_throw_MAPS_exception(errId);
elseif size(runDatasets,1) ~= size(provRunDatasets,1)
    errId = ['MAPS:',mfilename,':InconsistentOutputDims'];
    generate_and_throw_MAPS_exception(errId);
end

%% SETUP A MASTER EXCEPTION TO ADD CAUSES TO BELOW
% Setup a master exception to add causes to as encountered below.
masterErrId = ['MAPS:',mfilename,':BadCustomMacroOutputs'];
BadCustomMacroOutputsE = generate_MAPS_exception(masterErrId);

%% VALIDATE FORECAST & PROVISIONAL RUNS
% For each run output, validate the forecast and provisional run datasets
% output. Validate also that if the forecast run dataset output includes
% "past" data, then the provisional run dataset does also. This is so that
% EASE gets the right information to flag the runs output as "introduce
% news" as appropriate. If any of these validation checks fail, add the
% cause(s) of the failure(s) to the master exception.
errId = [masterErrId,':Instance'];
PreviousRunData = BaseRunData;
nRuns = size(runDatasets,1);
for iRun = 1:nRuns
    RunE = generate_MAPS_exception(errId,{num2str(iRun)});
    try
        validate_forecast_run_dataset_output(Model,runDatasets{iRun},T,H);
    catch RunContentValidationE
        RunE = addCause(RunE,RunContentValidationE);
    end
    try
        validate_provisional_run_dataset_output(...
            Model,provRunDatasets{iRun},T,H);
    catch ProvRunContentValidationE
        RunE = addCause(RunE,ProvRunContentValidationE);
    end
    RunE = validate_consistency_of_run_dataset_outputs(RunE,...
        Model,runDatasets{iRun},provRunDatasets{iRun},PreviousRunData);
    if ~isempty(RunE.cause)
        BadCustomMacroOutputsE = addCause(BadCustomMacroOutputsE,RunE);
    end
    PreviousRunData = runDatasets{iRun};
end

%% THROW ANY EXCEPTIONS ENCOUNTERED
% If any exceptions were encountered above the master exception will
% contain causes, so throw it back out into the EASE cutom macro interface.
if ~isempty(BadCustomMacroOutputsE.cause)
    throw(BadCustomMacroOutputsE);
end

end

%% FUNCTION TO VALIDATE THE CONISISTENCY OF THE RUN DATASETS
function RunE = validate_consistency_of_run_dataset_outputs(...
    RunE,Model,RunData,ProvRunData,PreviousRunData)
% This helper validates the past data outputs of a custom macro.
% It validates that the past data outputs in the provisional and forecast
% run dataset outputs are consistent with each other and the MAPS custom
% macro output rules.
%
% INPUTS:
%   -> RunE: Exception to add causes to
%   -> Model: MAPS model structure
%   -> RunData: forecast run dataset
%   -> ProvRunData: provisional run dataset
%   -> PreviousRunData: forecast run dataset from the previous run
%
% OUTPUTS:
%   -> RunE: Updated exception
%
% CALLS:
%   -> generate_MAPS_exception_and_add_as_cause
%   -> unpack_model
%   -> get_LSS_model_forecast_run_dataset_configs
%   -> get_NLBL_model_forecast_run_dataset_configs
%   -> cell_array_vlookup
%   -> unpack_data_from_MAPS_dataset

%% PAST DATA IS IN THE PROVISIONAL RUN BUT NOT THE RUN OUTPUT
% If past data is in the provisional run (implying that the run was an
% introduce news run) but not the forecast run, throw an exception.
if isfield(ProvRunData,'Past') && ~isfield(RunData,'Past')
    errId = ['MAPS:',mfilename,':PastDataInProvRunButNotInRun'];
    RunE = generate_MAPS_exception_and_add_as_cause(RunE,errId);
end

%% NEW PAST DATA IS IN THE RUN OUTPUT BUT NOT THE PROVISIONAL RUN
% If past data is in the forecast run output, but not the provisional run
% data and that past data is new (either because there was none before or
% it has changed), throw an exception.

if isfield(RunData,'Past') && ~isfield(ProvRunData,'Past')
    if ~isfield(PreviousRunData,'Past') || ...
            (isfield(PreviousRunData,'Past')&&...
            ~isequal(PreviousRunData.Past,RunData.Past))
        errId = ['MAPS:',mfilename,':NewPastDataInRunButNotInProvRun'];
        RunE = generate_MAPS_exception_and_add_as_cause(RunE,errId);
    end
end

%% PAST DATA IS DIFFERENT IN THE RUN OUTPUT & PROVISIONAL RUNS
% If past data exists in both the run & provisional run datasets, check
% that all overlapping datasets are identical. If any are not, throw an
% exception.
if isfield(RunData,'Past') && isfield(ProvRunData,'Past')
    modelIsLinearStateSpace = unpack_model(...
        Model,{'modelIsLinearStateSpace'});
    if modelIsLinearStateSpace
        runDataConfig = get_LSS_model_forecast_run_dataset_configs;
    else
        runDataConfig = get_NLBL_model_forecast_run_dataset_configs;
    end
    pastRunDataConfig = runDataConfig(strcmp(runDataConfig(:,1),'Past'),:);
    nPastRunDatasets = size(pastRunDataConfig,1);
    for iDataset = 1:nPastRunDatasets
        runPastDataset = unpack_data_from_MAPS_dataset(...
            RunData,pastRunDataConfig(iDataset,1:3));
        provRunPastDataset = unpack_data_from_MAPS_dataset(...
            ProvRunData,pastRunDataConfig(iDataset,1:3));
        if ~isempty(runPastDataset) && ~isempty(provRunPastDataset)
            if ~isequal(runPastDataset,provRunPastDataset)
                errId = ['MAPS:',mfilename,':InconsistentPastData'];
                RunE = generate_MAPS_exception_and_add_as_cause(...
                    RunE,errId);
                break
            end
        end
    end
end

%% CONSTRAINT DATA IS NOT CONSISTENT WITH PAST DATA IN THE RUN OUTPUT
% If past data and constraint data exists in the run dataset, check that
% they are consistent. If any are not, throw an exception.
if isfield(RunData,'Past') && isfield(RunData,'Constraint')
    modelIsLinearStateSpace = unpack_model(...
        Model,{'modelIsLinearStateSpace'});
    if modelIsLinearStateSpace
        runDataConfig = get_LSS_model_forecast_run_dataset_configs;
    else
        runDataConfig = get_NLBL_model_forecast_run_dataset_configs;
    end
    pastRunDataConfig = runDataConfig(strcmp(runDataConfig(:,1),'Past'),:);
    nPastRunDatasets = size(pastRunDataConfig,1);
    for iDataset = 1:nPastRunDatasets
        runPastDataset = unpack_data_from_MAPS_dataset(...
            RunData,pastRunDataConfig(iDataset,1:3));
        runConstraintDataset = unpack_data_from_MAPS_dataset(...
            RunData,[{'Constraint'} pastRunDataConfig(iDataset,2:3)]);
        if ~isempty(runPastDataset) && ~isempty(runConstraintDataset)
            nConstraintPeriods = size(runConstraintDataset,2);
            if any(any(abs(...
                    runPastDataset(:,end-nConstraintPeriods+1:end)-...
                    runConstraintDataset)>1e-12))
                errId = ['MAPS:',mfilename,...
                    ':InconsistentPastConstraintData'];
                RunE = generate_MAPS_exception_and_add_as_cause(...
                    RunE,errId);
                break
            end
        end
    end
end

end