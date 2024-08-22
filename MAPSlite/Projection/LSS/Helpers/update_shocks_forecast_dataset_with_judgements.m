function ShocksNew = update_shocks_forecast_dataset_with_judgements(...
    Model,Shocks,ShocksJudgements)
% This helper updates a LSS model shocks forecast dataset for judgements.
% It takes the values in the judgements input and overlays them on to the
% existing linear state space (LSS) model shocks forecast dataset.
%
% INPUTS:   
%   -> Model: MAPS LSS model structure
%   -> Shocks: structure of forecast shocks
%   -> ShocksJudgements: structure of judgements to forecast shock values
%
% OUTPUTS:  
%   -> ShocksNew: updated shocks forecast dataset
%
% CALLS:    
%   -> generate_and_throw_MAPS_exception
%   -> unpack_model
%   -> update_forecast_dataset_with_judgements
%
% DETAILS:  
%   -> This helper updates a linear state space model shocks forecast 
%      dataset for judgements.
%   -> It updates the forecast data values for each type of shock over
%      which judgement is being applied.
%   -> This function makes the assumption that the input shock judgements
%      and shocks are consistent with each other. If they are not, then it
%      will throw an exception but the cause of the failure may not be
%      captured precisely.
%           
% NOTES:    
%   -> Please see <> for a description of MAPS modules.
%   -> See also <> for a description of imposing judgement in MAPS.
%           
% This version: 04/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number and shape of the inputs is as expected. All inputs
% are compulsory.
if nargin < 3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);    
elseif ~isstruct(Shocks)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(ShocksJudgements)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
end
    
%% UNPACK SHOCKS MEATDATA FOR JUDGEMENTS
% Unpack the model meatdata for the shocks that is used in the judgement
% structure to specify which shock(s) the judgement(s) is to be applied to.
zMnems = unpack_model(Model,{'zMnems'});

%% SETUP OUTPUT
% Set the output as equal to the input to be updated below.
ShocksNew = Shocks;

%% DEFINE TYPES OF SHOCKS
% Define the types of shock that exist (anticipated & unanticipated) and 
% over which judgement could be applied. Compute the number of types to be
% used below.
anticipationTypes = {'anticipated';'unanticipated'};
nAnticipationTypes = size(anticipationTypes,1);

%% UPDATE SHOCKS
% Update each shock type over which judgement is being made. Catch any
% exceptions and add an appropriate message detailing which type of
% judgement could not be applied.
for iType = 1:nAnticipationTypes
    iAnticipationType = anticipationTypes{iType};
    if isfield(ShocksJudgements,iAnticipationType)
        try
            ShocksNew.(iAnticipationType) = ...
                update_forecast_dataset_with_judgements(zMnems,...
                Shocks.(iAnticipationType),...
                ShocksJudgements.(iAnticipationType));
        catch DatasetUpdateE
            errId = ['MAPS:',mfilename,':ShockDatasetUpdateFailure'];
            generate_MAPS_exception_add_cause_and_throw(...
                DatasetUpdateE,errId,{iAnticipationType});
        end
    end
end

end