function InversionInstructs = create_LSS_model_inversion_instructions(...
    Model,RunData,JudgementsInclude,JudgementsInfo)
% This LSS model macro helper creates the inputs required for inversion.
% It converts the fixes information in a linear state space (LSS) judgement
% information structure (the output of validation of a LSS model judgements
% dataset) into the format required for inversion of a LSS model to find
% values for the shocks that support the fixes.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> RunData: existing forecast run with (among other things):
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
%   -> JudgementsInclude: structure describing the judgements content
%       - anticipatedModelVariableFixes: true or false
%       - anticipatedModelObservableFixes: true or false
%       - anticipatedRawObservableFixes: true or false
%       - shocksForAnticipatedFixes: true or false
%       - unanticipatedModelVariableFixes: true or false
%       - unanticipatedModelObservableFixes: true or false
%       - unanticipatedRawObservableFixes: true or false
%       - shocksForUnanticipatedFixes: true or false
%       - anticipatedShockJudgements: true or false
%       - unanticipatedShockJudgements true or false
%       - timeVaryingTrendJudgements: true or false
%   -> JudgementsInfo: new structure with as many fields as indicated true
%      in the judgement content descriptor input, each of which contain:
%       - metadata: judgement metadata (from judgements dataset)
%       - cell: judgement data as a cell array (from judgements dataset)
%       - description: description of that type of judgement
%       - modelMetadata: complete model metadata for that variable type
%       - indices: model index numbers equivalent to the metadata
%       - matrix: matrix equivalent to the cell vectors
%
% OUTPUTS:  
%   -> InversionInstructs: structure describing new judgemental forecast
%       - AnticipatedFixes (optional): row cell array with columns denoting
%         periods (e.g. AnticipatedFixes.modelVarIndices{1} corresponds to
%         period 1 etc):
%           - modelVarIndices: indices of the model variables to be fixed
%           - modelVarValues: values for the fixes
%           - shockUsageIndices: indices of shocks to implement the
%             anticipated fixes
%       - UnanticipatedFixes (optional): as for anticipated fixes:
%           - modelVarIndices: indices of the model variables to be fixed
%           - modelVarValues: values for the fixes
%           - shockUsageIndices: indices of shocks to implement the
%             unanticipated fixes
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> convert_any_observable_fixes_to_model_variable_space (sub-function)
%   -> combine_all_fixes_in_model_space_together (sub-function)
%   -> convert_fix_info_to_inversion_instructions (sub-function)
%
% DETAILS:  
%   -> This linear state space (LSS) model forecast macro helper manages 
%      the translation of the judgement info output from judgement dataset 
%      validation to the format required for LSS model inversion.
%   -> It should be run after running the validation routine since it
%      consumes the output of that routine and because it assumes that the
%      inputs to this function are "correct" in the sense that they
%      constitute a valid set of judgements.
%   -> This helper first transforms all observable fixes (both raw and
%      model) to model variable space. 
%   -> It then combines all the fixes in model variable space together.
%   -> Finally, it translates them so that they are in the format required
%      for LSS model inversion.
%
% NOTES:
%   -> See the MAPS technical documentation for more details on this helper
%      and imposing judgement in LSS models.
%
% This version: 07/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that four inputs were passed in and then use a helper function to
% check that the incomplete model datasets and metadata passed in are
% valid.
if nargin < 4
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~isstruct(RunData)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId); 
elseif ~isstruct(JudgementsInclude)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(JudgementsInfo)
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);
end

%% CHECK JUDGEMENTS INPUT IS CONSISTENT WITH THIS FUNCTION
% The judgement info structure must include fix information or this
% function cannot produce an inversion instructions structure.
if ~isfield(JudgementsInclude,'anticipatedShockJudgements') && ...
        ~isfield(JudgementsInclude,'unanticipatedShockJudgements')
    errId = ['MAPS:',mfilename,':BadJudgementsInput'];
    generate_and_throw_MAPS_exception(errId);    
end

%% CONVERT OBSERVABLE FIXES TO MODEL VARIABLES SPACE
% Call a sub-function to convert any fixes made in observable space into
% model variable space. This function returns a structure with the same
% fields as the judgements information structure input, but with the
% content converted to model variable space.
JudgementsInfoInModelSpace = ...
    convert_any_observable_fixes_to_model_variable_space(...
    Model,RunData,JudgementsInclude,JudgementsInfo);

%% COMBINE ALL ANTICIPATED & UNANTICIPATED FIXES IN MODEL SPACE TOGETHER
% Call a sub-function to combine all fixes made under a particular 
% assumption together in model variable space.
CombinedFixes = combine_all_fixes_in_model_space_together(...
    JudgementsInclude,JudgementsInfoInModelSpace);

%% CREATE INVERSION INSTRUCTIONS
% Call another sub-function to create the inversion instructions from the 
% fixes information in the judgement structure. 
InversionInstructs = convert_fix_info_to_inversion_instructions(...
    CombinedFixes);

end

%% FUNCTION TO CREATE INVSERSION INSTRUCTIONS
function InversionInstructs = ...
    convert_fix_info_to_inversion_instructions(CombinedFixes)
% This helper function creates all the inputs required for inversion.
% It transforms the judgements information into the format required for the 
% inversion module. The main transformation (in addition to getting the 
% structure of the input right) is for fixes to be organised by time 
% periods (instead of by variable as they were in the input to this 
% function).
%
% INPUTS:   
%   -> CombinedFixes: structure of fix judgement information
%
% OUTPUTS:  
%   -> InversionInstructs: structure of information required for inversion
%
% CALLS:
%   -> create_time_varying_cell_arrays (sub-function)

%% RUN THROUGH ANTICIPATED % UNANTICIPATED FIXES
% Run through the two fix types in turn. If they exist in the fixes input,
% then convert both the indices and matrices in the fixes and shocks to
% cell arrays organised by time rather than by variable. Once the shocks
% have been converted, reflect the shock set so that the indices and values
% refer to shocks not being used rather than the shocks being used.
fixTypes = {'Anticipated';'Unanticipated'};
nFixTypes = size(fixTypes,1);
for iType = 1:nFixTypes
    iFixType = fixTypes{iType};
    iFixTypeName = [iFixType,'Fixes'];
    if isfield(CombinedFixes,iFixType)
        [InversionInstructs.(iFixTypeName).modelVarIndices,...
            InversionInstructs.(iFixTypeName).modelVarValues] = ...
            create_time_varying_cell_arrays(...
            CombinedFixes.(iFixType).ModelVar.indices,...
            CombinedFixes.(iFixType).ModelVar.matrix);
        InversionInstructs.(iFixTypeName).shockUsageIndices = ...
            create_time_varying_cell_arrays(...
            CombinedFixes.(iFixType).Shocks.indices,...
            CombinedFixes.(iFixType).Shocks.matrix);
    end
end

end

%% FUNCTION TO CREATE CELL ARRAY OF NAMES & VALUES
function [jCellInd,jCell] = create_time_varying_cell_arrays(jInd,jMat)
% This helper function converts fix matrix information into cell arrays.
% It forms part of inversion input creation as it creates judgement
% cell array inputs organised by time periods as required in the inversion 
% algorithm.
%
% INPUTS:   
%   -> jInd: vector of model index numbers
%   -> jMat: numerical matrix of judgements of dimension size(jInd,1)*S
%
% OUTPUTS:  
%   -> jCellInd: 1*S cell array of model index vectors
%   -> jCell: 1*S cell array of numerical judgemnt vectors
%
% CALLS:
%   -> none

%% SETUP OUTPUT
% Both the index and value outputs are a cell array of vectors of the same
% length as number of inversion periods.
S = size(jMat,2);
jCellInd = cell(1,S);
jCell = cell(1,S);

%% COMPUTE THE NON-NAN LOGICALS IN THE JUDGEMENT MATRIX
% Compute a matrix of logicals which indicate the elements of the judgement
% matrix that are non-nan.
notNaNlogicals = ~isnan(jMat);

%% COMPUTE THE SHOCK INDICES & VALUES PERIOD-BY PERIOD
% Run-through the inversion periods setting the index numbers and logicals
% in accordance with the non-NaN values computed above.
for s = 1:S
    jCellInd{s} = jInd(notNaNlogicals(:,s));
    jCellInd{s} = jCellInd{s}(:);
    jCell{s} = jMat(notNaNlogicals(:,s),s);
end

end

%% FUNCTION TO COMBINE ALL FIXES IN MODEL SPACE TOGETHER
function CombinedFixes = combine_all_fixes_in_model_space_together(...
    JudgementsInclude,JudgementsInfoInModelSpace)
% This helper combines all fixes together in model variable space.
% It combines all of the fixes made under each assumption together ready
% for conversion to the format required for inversion (eg it combines any
% anticipated raw observable fixes and model observable fixes already
% converted to model variable space with any anticipated model variable
% fixes included in the input.
%
% INPUTS:   
%   -> JudgementsInclude: structure describing the judgements content
%       - anticipatedModelVariableFixes: true or false
%       - anticipatedModelObservableFixes: true or false
%       - anticipatedRawObservableFixes: true or false
%       - shocksForAnticipatedFixes: true or false
%       - unanticipatedModelVariableFixes: true or false
%       - unanticipatedModelObservableFixes: true or false
%       - unanticipatedRawObservableFixes: true or false
%       - shocksForUnanticipatedFixes: true or false
%       - anticipatedShockJudgements: true or false
%       - unanticipatedShockJudgements true or false
%       - timeVaryingTrendJudgements: true or false
%   -> JudgementsInfoInModelSpace: structure with as many fields as 
%      indicated true in the judgement content descriptor input, each of 
%      which contain (among other things):
%       - indices: model index numbers in model variable space
%       - matrix: potentially incomplete matrix of fix data
%
% OUTPUTS:  
%   -> CombinedFixes: combined fixes in model variable space:
%       - Anticipated: anticipated fixes:
%           - ModelVar: fixes on the model variables:
%               - indices: model index numbers for these fixes
%               - matrix: matrix of values for the fixes
%           - Shocks: shocks to be used to achieve the fixes:
%               - indices: model index numbers for these shocks
%               - matrix: matrix of usage indicators for the shocks
%       - Unanticipated: anticipated fixes:
%           - ModelVar: fixes on the model variables:
%               - indices: model index numbers for these fixes
%               - matrix: matrix of values for the fixes
%           - Shocks: shocks to be used to achieve the fixes:
%               - indices: model index numbers for these shocks
%               - matrix: matrix of usage indicators for the shocks
%
% CALLS:
%   -> overlay_one_incomplete_dataset_on_to_another

%% SETUP OUTPUT
% Setup the output to be updated below.
CombinedFixes = struct;

%% COMBINE FIXES
% Define the fix and space types. For each combination, check whether or
% not that fix type exists in the judgements input. If it does, either 
% overlay it on to an existing fix of that type (i.e. combine fixes 
% together) or create a new one.  
fixTypes = {'Anticipated';'Unanticipated'};
nFixTypes = size(fixTypes,1);
for iType = 1:nFixTypes
    spaceTypes = {'RawObservable';'ModelObservable';'ModelVariable'};
    nSpaceTypes = size(spaceTypes,1);
    for iiType = 1:nSpaceTypes
        iFixSpaceType = ...
            [lower(fixTypes{iType}),spaceTypes{iiType},'Fixes'];
        if JudgementsInclude.(iFixSpaceType)
            if isfield(CombinedFixes,fixTypes{iType})
                [CombinedFixes.(fixTypes{iType}).ModelVar.matrix,...
                    CombinedFixes.(fixTypes{iType}).ModelVar.indices] = ...
                    overlay_one_incomplete_dataset_on_to_another(...
                    CombinedFixes.(fixTypes{iType}).ModelVar.matrix,...
                    CombinedFixes.(fixTypes{iType}).ModelVar.indices,...
                    JudgementsInfoInModelSpace.(iFixSpaceType).matrix,...
                    JudgementsInfoInModelSpace.(iFixSpaceType).indices);               
            else
                CombinedFixes.(fixTypes{iType}).ModelVar.matrix = ...
                    JudgementsInfoInModelSpace.(iFixSpaceType).matrix;
                CombinedFixes.(fixTypes{iType}).ModelVar.indices = ...
                    JudgementsInfoInModelSpace.(iFixSpaceType).indices;                
            end
        end
    end
    
end

%% ADD SHOCK USAGE DATA
% Add the shock usage data to the combined fixes structure as appropriate
% given their existence or otherwise in the input.
shockTypes = {'Anticipated';'Unanticipated'};
nShockTypes = size(shockTypes,1);
for iType = 1:nShockTypes
    iShockType = ['shocksFor',fixTypes{iType},'Fixes'];
    if JudgementsInclude.(iShockType)        
        CombinedFixes.(shockTypes{iType}).Shocks.matrix = ...
            JudgementsInfoInModelSpace.(iShockType).matrix;
        CombinedFixes.(shockTypes{iType}).Shocks.indices = ...
            JudgementsInfoInModelSpace.(iShockType).indices;
        
    end 
end

end

%% FUNCTION TO CONVERT OBSERVABLE FIXES TO MODEL VARIABLE SPACE
function JudgementsInfoInModelSpace = ...
    convert_any_observable_fixes_to_model_variable_space(Model,...
    RunData,JudgementsInclude,JudgementsInfo)
% This helper converts any obseravble fixes to model variable space.
% It uses the data transformation symbolics and measurement equation 
% matrices in two sub-functions to convert any observable fixes included in 
% the judgement info structure to model variable space.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> RunData: existing forecast run with (among other things):
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
%   -> JudgementsInclude: structure describing the judgements content
%       - anticipatedModelVariableFixes: true or false
%       - anticipatedModelObservableFixes: true or false
%       - anticipatedRawObservableFixes: true or false
%       - shocksForAnticipatedFixes: true or false
%       - unanticipatedModelVariableFixes: true or false
%       - unanticipatedModelObservableFixes: true or false
%       - unanticipatedRawObservableFixes: true or false
%       - shocksForUnanticipatedFixes: true or false
%       - anticipatedShockJudgements: true or false
%       - unanticipatedShockJudgements true or false
%       - timeVaryingTrendJudgements: true or false
%   -> JudgementsInfo: structure with as many fields as 
%      indicated true in the judgement content descriptor input, each of 
%      which contain (among other things):
%       - indices: model index numbers in model variable space
%       - matrix: potentially incomplete matrix of fix data
%
% OUTPUTS:  
%   -> JudgementsInfoInModelSpace: updated version of the input with all
%      observable fixes converted into model variable space
%       - indices: model index numbers in model variable space
%       - matrix: potentially incomplete matrix of fix data
%
% CALLS:
%   -> convert_both_raw_observable_fixes_to_model_observable_space
%      (sub-function)
%   -> convert_raw_observable_fixes_to_model_observable_space 
%      (sub-function)
%   -> unpack_model
%   -> find_index_numbers_of_observables_in_model_var_space
%   -> convert_model_observable_fixes_to_model_variable_space
%      (sub-function)

%% SETUP OUTPUT
% Setup the output to this function as equal to the input. This structure
% will be updated in the rest of this function if the judgements contain
% any observable fixes.
JudgementsInfoInModelSpace = JudgementsInfo;

%% CONVERT RAW OBSERVABLE FIXES TO MODEL OBSERAVABLE SPACE IF BOTH EXIST
% If both anticipated and unanticipated fixes exist in the judgements info,
% then the conversion is slightly different. It is first necessary to
% combine the judgements together before converting them to model
% observable space. One of the rules of impose judgement (associated with
% lags that may exist in the data transformation equations) is that
% judgements made in raw observable space must be continuous. It is
% possible to satisfy that rule through a combination of anticipated and
% unanticipated judgements, so it is first necessary to combine those
% together before transforming the data to guarantee that the result is
% valid (given that the combined set of judgements meets the rules).
% And, because it is necessary to keep track of anticipated and
% unanticipated judgements seprately, it is then necessary to transform the
% result afterwards.
if JudgementsInclude.anticipatedRawObservableFixes && ...
        JudgementsInclude.unanticipatedRawObservableFixes
    fixType1 = 'anticipatedRawObservableFixes';
    fixType2 = 'unanticipatedRawObservableFixes';
    [JudgementsInfoInModelSpace.(fixType1).matrix,...
        JudgementsInfoInModelSpace.(fixType1).indices,...
        JudgementsInfoInModelSpace.(fixType2).matrix,...
        JudgementsInfoInModelSpace.(fixType2).indices] = ...
        convert_both_raw_observable_fixes_to_model_observable_space(...
        Model,RunData,...
        JudgementsInfo.(fixType1).matrix,...
        JudgementsInfo.(fixType1).indices,...
        JudgementsInfo.(fixType2).matrix,...
        JudgementsInfo.(fixType2).indices);
end
    
%% CONVERT RAW OBSERVABLE FIXES TO MODEL OBSERAVABLE SPACE IF ONE EXISTS
% If the judgements only include anticipated or unanticipated raw 
% observable fixes (or neither) then they can be transformed to model 
% observable space individually.
if (JudgementsInclude.anticipatedRawObservableFixes&&...
        ~JudgementsInclude.unanticipatedRawObservableFixes) || ...
        (~JudgementsInclude.anticipatedRawObservableFixes&&...
        JudgementsInclude.unanticipatedRawObservableFixes)
    fixTypes = {'anticipated';'unanticipated'};
    nFixTypes = size(fixTypes,1);
    for iType = 1:nFixTypes
        iFixType = [fixTypes{iType},'RawObservableFixes'];
        if JudgementsInclude.(iFixType)
            [JudgementsInfoInModelSpace.(iFixType).matrix,...
                JudgementsInfoInModelSpace.(iFixType).indices] = ...
                convert_raw_observable_fixes_to_model_observable_space(...
                Model,RunData,JudgementsInfo.(iFixType).matrix,...
                JudgementsInfo.(iFixType).indices);
        end
    end
end

%% CONVERT ALL FIXES IN MODEL OBSERVABLE SPACE TO MODEL VARIABLE SPACE
% Convert all raw observable fixes that were converted to model observable 
% space above and all model observable fixes to model variable space using
% the helper function below.
if JudgementsInclude.anticipatedRawObservableFixes || ... 
    JudgementsInclude.unanticipatedRawObservableFixes || ...
    JudgementsInclude.anticipatedModelObservableFixes || ...
    JudgementsInclude.unanticipatedModelObservableFixes
    YxInds = find_index_numbers_of_observables_in_model_var_space(Model);
    [D,G] = unpack_model(Model,{'D','G'});
    fixTypes = {'anticipatedRawObservable'
        'unanticipatedRawObservable'
        'anticipatedModelObservable'
        'unanticipatedModelObservable'};
    nFixTypes = size(fixTypes,1);
    for iType = 1:nFixTypes
        iFixType = [fixTypes{iType},'Fixes'];            
        if JudgementsInclude.(iFixType)
            [JudgementsInfoInModelSpace.(iFixType).matrix,...
                JudgementsInfoInModelSpace.(iFixType).indices] = ...
                convert_model_observable_fixes_to_model_variable_space(...
                D,G,YxInds,...
                JudgementsInfoInModelSpace.(iFixType).matrix,...
                JudgementsInfoInModelSpace.(iFixType).indices);
        end
    end
end

end

%% FUNCTION TO COMBINE RAW FIXES & CONVERT TO MODEL OBSERVABLE SPACE
function [YtildeYaMat,YtildeYaInds,YtildeYuMat,YtildeYuInds] = ...
    convert_both_raw_observable_fixes_to_model_observable_space(...
    Model,RunData,YtildeaMat,YtildeaInds,YtildeuMat,YtildeuInds)
% This helper converts both sets of raw obseravble fixes to model space.
% It combines the raw observable fixes being made under both the
% anticipated and unanticipated assumptions together, converts the combined
% set of raw observable fixes to model space and then separates them back
% out into transformed anticipated and unanticipated fixes.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> RunData: existing forecast run with (among other things):
%       - Constraint: structure of forecast initial conditions
%           - rawObservables: vector of raw observable 
%             constraints
%           - timeVaryingTrends (model dep.): vector of time-varying trend 
%             constraints
%       - Forecast: stucture of forecast data with (among other things):
%           - timeVaryingTrends (model dep.): matrix of time-varying trends
%             over the forecast
%   -> YtildeaMat: incomplete anticipated raw observable judgements dataset
%   -> YtildeaInds: model index numbers for the anticipated raw observables 
%      being fixed
%   -> YtildeuMat: incomplete unanticipated raw observable judgements 
%      dataset
%   -> YtildeuInds: model index numbers for the unanticipated raw 
%      observables being fixed
%
% OUTPUTS:  
%   -> YtildeYaMat: incomplete anticipated raw observable judgements 
%      dataset transformed to model observable space
%   -> YtildeYaInds: model index numbers for the transformed anticipated 
%      raw observable judgements
%   -> YtildeYuMat: incomplete unanticipated raw observable judgements 
%      dataset transformed to model observable space
%   -> YtildeYuInds: model index numbers for the transformed unanticipated 
%      raw observable judgements
%
% CALLS:
%   -> overlay_one_incomplete_dataset_on_to_another
%   -> convert_raw_observable_fixes_to_model_observable_space
%   (sub-function)
%   -> undo_overlay_of_one_incomplete_dataset_on_to_another

%% COMBINE ANTICIPATED & UNANTICIPATED FIXES TOGETHER
% Use a MAPS helper function to combine the anticipated & unanticipated raw
% observable fixes together. This is so that the raw observable fixes can 
% be transformed into model observable space jointly to allow for the fact
% that a valid, continuous set of raw observable fixes can be made through
% a combination fixes made under the anticipated & unanticipated 
% assumptions. See the MAPS rules on imposing judgement in linear state
% space (LSS) models for more details.
[YtildejMat,YtildejInds,YtildeaLogicals,YtildeuLogicals] = ...
    overlay_one_incomplete_dataset_on_to_another(...
    YtildeaMat,YtildeaInds,YtildeuMat,YtildeuInds);

%% COMBINE COMBINED RAW OBSERVABLE FIXES TO MODEL OBSERVABLE SPACE
% Use the sub-function below to transform the combined set of fixes into
% model observable space.
[YtildeYjMat,YtildeYjInds] = ...
    convert_raw_observable_fixes_to_model_observable_space(...
    Model,RunData,YtildejMat,YtildejInds);

%% SEPARATE THE CONVERTED OBSERVABLE FIXES INTO ANTICIPATED & UNANTICIPATED
% Separate the transformed set of fixes back out into its constituent 
% anticipated & unanticipated parts. 
[YtildeYaMat,YtildeYaInds,YtildeYuMat,YtildeYuInds] = ...
    undo_overlay_of_one_incomplete_dataset_on_to_another(...
    YtildeYjMat,YtildeYjInds,YtildeaLogicals,YtildeuLogicals);

end

%% FUNCTION TO CONVERT RAW OBSERVABLE FIXES TO MODEL OBSERVABLE SPACE
function [YtildeYjMat,YtildeYjInds] = ...
    convert_raw_observable_fixes_to_model_observable_space(...
    Model,RunData,YtildejMat,YtildejInds)
% This helper converts arbitrary raw obseravble fixes to model space.
% It uses the data transformation symbolics to convert any observable 
% fixes included in the judgement info structure to model observable space.
% It then validates the result.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> RunData: existing forecast run with (among other things):
%       - Constraint: structure of forecast initial conditions
%           - rawObservables: vector of raw observable 
%             constraints
%           - timeVaryingTrends (model dep.): vector of time-varying trend 
%             constraints
%       - Forecast: stucture of forecast data with (among other things):
%           - timeVaryingTrends (model dep.): matrix of time-varying trends
%             over the forecast
%   -> YtildejMat: incomplete raw observable judgements dataset
%   -> YtildejInds: model index numbers for the variables being fixed
%
% OUTPUTS:  
%   -> YtildeYjMat: incomplete raw observable judgements dataset
%      transformed to model observable space
%   -> YtildeYjInds: model index numbers for the transformed raw observable
%      judgements
%
% CALLS:
%   -> unpack_model
%   -> transform_indices_from_raw_to_model_observable_space
%   -> transform_observables_from_raw_to_model_space
%   -> validate_transformed_model_observable_judgements

%% UNPACK RUN DATA & MODEL
% Unpack the data transformation function handle, the time trend indicator
% and model observable mnemonics from the model input. Unpack the raw
% observable forecast initial conditions and, if the model has time trends,
% the time-varying trend initial conditions and forecast data. Note that
% this cell assumes that the combination of judgements, data and model is
% valid (i.e. that they have been through the judgements dataset validation
% function).
[DTfunHandle,modelHasTimeVaryingTrends,Ymnems] = unpack_model(...
    Model,{'DTfunHandle','modelHasTimeVaryingTrends','Ymnems'});
YtildeT = RunData.Constraint.rawObservables;
if modelHasTimeVaryingTrends
    etatT = RunData.Constraint.timeVaryingTrends;
    etatf = RunData.Forecast.timeVaryingTrends;
end

%% CREATE A TEMPORARY FULL MATRIX OF RAW OBSERVABLES
% Create a temporary full matrix of raw observables with NaNs inserted for
% the observables not being fixed.
SYj = size(YtildejMat,2);
nY = size(YtildeT,1);
YtildejMatTemp = NaN*ones(nY,SYj);
YtildejMatTemp(YtildejInds,:) = YtildejMat;

%% COMPUTE INDEX NUMBERS OF RAW OBSERVABLE FIXES IN MODEL OBSERVABLE SPACE
% Call a helper function to compute the index numbers of the raw
% observables in model observable space. ompute the metedata associated
% with those index numbers for validation below.
YtildeYjInds = transform_indices_from_raw_to_model_observable_space(...
    Model,YtildejInds);
YtildeYjMetadata = Ymnems(YtildeYjInds);

%% TRANSFORM THE DATA
% Use the MAPS data transformation function to convert the fixes in raw
% observable space into model observable space. Note that the function 
% retains the NaN indicator for no usage in a particular period. However,
% the validity of the result relies on the raw observable fixes meeting the
% impose judgement rules (because data transformation equations with lags
% in them will carry lagged NaN values through into observable space so raw
% observable fixes must be continuous).
if modelHasTimeVaryingTrends
    YtildeYjMatTemp = transform_observables_from_raw_to_model_space(...
        DTfunHandle,YtildejMatTemp,YtildeT,etatf(:,1:SYj),etatT);
else
    YtildeYjMatTemp = transform_observables_from_raw_to_model_space(...
        DTfunHandle,YtildejMatTemp,YtildeT);
end
YtildeYjMat = YtildeYjMatTemp(YtildeYjInds,:);

%% VALIDATE THE TRANSFORMED DATA
% Find any periods at the end of the raw observable judgements dataset
% input to this function that were all NaNs (since these represent periods
% where no raw observable fixes are being applied). Validate the remaining
% transformed data to check that it does not have imaginary, inf, or NaN
% values. This function will throw an exception detailing the metadata
% associated with the series that fail validation. The purpose of this
% check is to validate the combination of the raw obseravable fix data, the
% time trends and the data transformation equations.
nYj = size(YtildeYjMat,1);
SYjs = SYj*ones(nYj,1);
[indYj,SiYjs] = find(diff(isnan(YtildejMat),1,2),nYj,'last');
SYjs(indYj) = SiYjs;
validate_transformed_model_observable_judgements(...
    YtildeYjMat,YtildeYjMetadata,SYjs);

end

%% FUNCTION TO CONVERT MODEL OBSERVABLE FIXES TO MODEL VARIABLE SPACE
function [YxjMat,YxjInds] = ...
    convert_model_observable_fixes_to_model_variable_space(...
    D,G,YxInds,YjMat,YjInds)
% This helper converts fixes in model observable space to model variable 
% space in preparation for model inversion.
% 
% INPUTS:   
%   -> D: constants in the measurement equation
%   -> G: loadings on the model variables in the measurement equations
%   -> YxInds: indices of model observables in model variable space
%   -> YjMat: model observables judgements matrix
%   -> YjInd: indexes of judgements in model observable space
%
% OUTPUTS:  
%   -> YxjMat: observable judgements data translated to model space
%   -> YxjInds: indexes of the observables in model variable space
%
% CALLS:
%   -> none

%% COMPUTE THE INDICES OF THE FIXES IN MODEL VARIABLE SPACE
% Use the input index of model observables in model variable space and the
% indices of the fixes in model observable space to compute the index
% numbers of the observable fixes transformed to model variable space.
YxjInds = YxInds(YjInds);

%% COMPUTE THE MATRIX OF JUDGEMENTS IN MODEL VARIABLE SPACE
% Invert the measurement equation to compute the matrix of judgements in
% model variable space.
Sj = size(YjMat,2);
YxjMat = (YjMat-D(YjInds)*ones(1,Sj))...
    ./(diag(G(YjInds,YxjInds))*ones(1,Sj));

end