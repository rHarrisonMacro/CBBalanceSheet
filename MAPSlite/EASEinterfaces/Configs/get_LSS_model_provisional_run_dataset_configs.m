function [MAPSprovRunDataConfig,EASEprovRunDataConfig] = ...
    get_LSS_model_provisional_run_dataset_configs
% This config contains MAPS/EASE LSS model prov. run dataset formats. 
% It contains information about the format of linear state space (LSS) 
% model provisional run datasets passed from EASE to MAPS and information 
% about how the information is represented in MAPS. Provisional run 
% information is information that a user wishes to update the projection 
% with (eg new back data or new judgements).
%
% INPUTS:   
%   -> none
%
% OUTPUTS: 
%   -> MAPSprovRunDataConfig: a cell array of config info for MAPS-
%      represented LSS model forecast run datasets and how they relate to 
%      the model
%   -> EASEprovRunDataConfig: a cell array of config info for EASE-
%      represented LSS model provisional run datasets
%
% CALLS:
%   -> none
%
% DETAILS:
%   -> Data passed into MAPS from EASE is represented in n-column cell
%      arrays. The first n-2 columns contain metadata about the type of
%      data; the penultimate column contains metadata about the specific 
%      series (which links it to a variable in the model); the final column
%      contains a row time series vector of data for the variable.
%   -> MAPS represents the data in a more natural matrix oriented modelling
%      environment, in which complete datasets for a particular variable
%      type are stored in matrices (of dimension nVar*T). Each of these
%      matrices in stored in a single structure whose fieldnames play the 
%      same role as the metadata in the first n-2 columns of the EASE
%      represented dataset.
%   -> For example, a complete dataset for the raw observables over the
%      past (i.e. back data) in a model with 20 observables would occupy 20 
%      rows of a five column cell array passed to MAPS from EASE. The first
%      column would contain the string 'Past' in each of the 20 elements;
%      the second would contain the string 'rawObservables' in each row;
%      the third would be empty (because those two strings are sufficient
%      to identify the data type); the fourth would contain model mnemonics
%      - one for each of the 20 observables; the fifth would contain a
%      vector of data points in each element for each of the variables. In
%      MAPS this data would be combined into a matrix and reordered to 
%      ensure consistency with the model order for raw observables and 
%      would be stored in a structure with field Past.rawObservables.
%   -> MAPS stores incomplete datasets (i.e. where the dataset does not 
%      include all the variables in the model, such as data describing
%      judgements) the same way (i.e. using a structure) but represents 
%      them in 2 column cell arrays with the first column containing model
%      metadata to associate the data with a particular variable in the
%      model and the second column containing individual data series for
%      the variables included (i.e. the same format as in the final two
%      columns of EASE-represented datasets).
%
% NOTES:
%   -> See xxxxxx for a description of EASE interfaces in MAPS and how the 
%      data translation is managed. 
%
% This version: 14/03/2011
% Author(s): Matt Waldron

%% DEFINE MAPS LSS MODEL PROVISIONAL RUN DATASET CONFIGURATION
% Define the information that describes how MAPS represents linear state 
% space model provisional run datasets (under the assumption that MAPS 
% datasets are stored as structures as described above). Define also how 
% the data relates to the model being used, whether the information is
% compulsory or optional and the number of data points that should exist in
% datasets for that data type (where H is the forecast horizon and T is the
% past data horizon in model variable space).
MAPSprovRunDataConfig = {
    'Forecast'      'UnanticipatedFixes'    'modelVariables'    'cell'      'numeric'               'xMnems'        'optional'                                      'H'
    'Forecast'      'AnticipatedFixes'      'modelVariables'    'cell'      'numeric'               'xMnems'        'optional'                                      'H'
    'Forecast'      'UnanticipatedFixes'    'modelObservables'  'cell'      'numeric'               'Ymnems'        'optional'                                      'H'
    'Forecast'      'AnticipatedFixes'      'modelObservables'  'cell'      'numeric'               'Ymnems'        'optional'                                      'H'
    'Forecast'      'UnanticipatedFixes'    'rawObservables'    'cell'      'numeric'               'YtildeMnems'   'optional'                                      'H'
    'Forecast'      'AnticipatedFixes'      'rawObservables'    'cell'      'numeric'               'YtildeMnems'   'optional'                                      'H'
    'Forecast'      'Shocks'                'unanticipated'     'cell'      'numeric'               'zMnems'        'optional'                                      'H'
    'Forecast'      'Shocks'                'anticipated'       'cell'      'numeric'               'zMnems'        'optional'                                      'H'
    'Forecast'      'timeVaryingTrends'     ''                  'cell'      'numeric'               'etatMnems'     'optional'                                      'H'
    'Forecast'      'UnanticipatedFixes'    'shockUsages'       'cell'      'logical equivalent'    'zMnems'        'optional'                                      'H'
    'Forecast'      'AnticipatedFixes'      'shockUsages'       'cell'      'logical equivalent'    'zMnems'        'optional'                                      'H'
    'Past'          'rawObservables'        ''                  'matrix'    'numeric'               'YtildeMnems'   'compulsory if modelHasDataTransformationEqs'   'T+1'
    'Past'          'timeVaryingTrends'     ''                  'matrix'    'numeric'               'etatMnems'     'compulsory if modelHasTimeVaryingTrends'       'T+1'
    };

%% DEFINE EASE LSS MODEL PROVISIONAL RUN DATASET CONFIGURATION
% Define the information that describes how EASE represents LSS model 
% provisional run datasets (under the assumption that EASE datasets are 
% stored as cell arrays as described in the notes above).
EASEprovRunDataConfig = {
    'Fixes'         'ModelVariables'        'unanticipated'
    'Fixes'         'ModelVariables'        'anticipated'
    'Fixes'         'ModelObservables'      'unanticipated'
    'Fixes'         'ModelObservables'      'anticipated'
    'Fixes'         'RawObservables'        'unanticipated'
    'Fixes'         'RawObservables'        'anticipated'
    'Fixes'         'Shocks'                'unanticipated'
    'Fixes'         'Shocks'                'anticipated'
    'Fixes'         'TimeVaryingTrends'     ''
    'Usages'        'Shocks'                'unanticipated'
    'Usages'        'Shocks'                'anticipated'
    'Past'          'RawObservables'        ''
    'Past'          'TimeVaryingTrends'     ''
    };

end