function [MAPSrunDataConfig,EASErunDataConfig] = ...
    get_LSS_model_forecast_run_dataset_configs
% This config contains MAPS/EASE LSS model forecast run dataset formats. 
% It contains information about the format of linear state space model 
% forecast run datasets passed from EASE to MAPS and information about how
% to translate that to a MAPS formatted equivalent.
%
% INPUTS:   
%   -> none
%
% OUTPUTS:  
%   -> MAPSrunDataConfig: a cell array of config info for MAPS-represented
%      LSS model forecast run datasets and how they relate to the model
%   -> EASErunDataConfig: a cell array of config info for EASE-represented
%      LSS model forecast run datasets
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

%% LIST MAPS LSS MODEL FORECAST RUN DATASET CONFIGURATION
% Define the information that describes how MAPS represents linear state 
% space model forecast run datasets (under the assumption that MAPS 
% datasets are stored as structures as described above). Define also how 
% the data relates to the model being used, whether the information is
% compulsory or optional and the number of data points that should exist in
% datasets for that data type (where H is the forecast horizon and T is the
% past data horizon in model variable space).
MAPSrunDataConfig = {
    'Past'          'modelVariables'        ''                  'matrix'    'numeric'   'xMnems'        'compulsory'                                    'T'
    'Past'          'Shocks'                'unanticipated'     'matrix'    'numeric'   'zMnems'        'compulsory'                                    'T'
    'Past'          'modelObservables'      ''                  'matrix'    'numeric'   'Ymnems'        'compulsory if modelHasMeasurementEqs'          'T'
    'Past'          'measurementErrors'     ''                  'matrix'    'numeric'   'wMnems'        'compulsory if modelHasMeasurementErrors'       'T'
    'Past'          'rawObservables'        ''                  'matrix'    'numeric'   'YtildeMnems'   'compulsory if modelHasDataTransformationEqs'   'T+1'
    'Past'          'timeVaryingTrends'     ''                  'matrix'    'numeric'   'etatMnems'     'compulsory if modelHasTimeVaryingTrends'       'T+1'
    'Constraint'    'modelVariables'        ''                  'matrix'    'numeric'   'xMnems'        'compulsory'                                    '1'
    'Constraint'    'rawObservables'        ''                  'matrix'    'numeric'   'YtildeMnems'   'compulsory if modelHasDataTransformationEqs'   '1'
    'Constraint'    'timeVaryingTrends'     ''                  'matrix'    'numeric'   'etatMnems'     'compulsory if modelHasTimeVaryingTrends'       '1'
    'Forecast'      'modelVariables'        ''                  'matrix'    'numeric'   'xMnems'        'compulsory'                                    'H'
    'Forecast'      'Shocks'                'unanticipated'     'matrix'    'numeric'   'zMnems'        'compulsory'                                    'H'
    'Forecast'      'Shocks'                'anticipated'       'matrix'    'numeric'   'zMnems'        'compulsory if modelIsForwardLooking'           'H'
    'Forecast'      'modelObservables'      ''                  'matrix'    'numeric'   'Ymnems'        'compulsory if modelHasMeasurementEqs'          'H'
    'Forecast'      'measurementErrors'     ''                  'matrix'    'numeric'   'wMnems'        'optional'                                      'H'
    'Forecast'      'rawObservables'        ''                  'matrix'    'numeric'   'YtildeMnems'   'compulsory if modelHasDataTransformationEqs'   'H'
    'Forecast'      'timeVaryingTrends'     ''                  'matrix'    'numeric'   'etatMnems'     'compulsory if modelHasTimeVaryingTrends'       'H'
    };

%% DEFINE EASE LSS MODEL FORECAST RUN DATASET CONFIGURATION
% Define the information that describes how EASE represents LSS model 
% forecast run datasets (under the assumption that EASE datasets are 
% stored as cell arrays as described in the notes above).
EASErunDataConfig = {
    'Past'          'modelVariables'        ''
    'Past'          'Shocks'                'unanticipated'
    'Past'          'modelObservables'      ''
    'Past'          'measurementErrors'     ''
    'Past'          'rawObservables'        ''
    'Past'          'timeVaryingTrends'     ''
    'Constraint'    'modelVariables'        ''
    'Constraint'    'rawObservables'        ''
    'Constraint'    'timeVaryingTrends'     ''
    'Forecast'      'modelVariables'        ''
    'Forecast'      'Shocks'                'unanticipated'
    'Forecast'      'Shocks'                'anticipated'
    'Forecast'      'modelObservables'      ''
    'Forecast'      'measurementErrors'     ''
    'Forecast'      'rawObservables'        ''
    'Forecast'      'timeVaryingTrends'     ''
    };

end