function [EASErunData,EASEprovRunData] = ...
    convert_run_datasets_from_MAPS_to_EASE_format(...
    Model,MAPSrunData,MAPSprovRunData)
% This helper converts forecast & prov run datasets from MAPS to EASE.
% It should be used as a helper to all EASE interfaces in MAPS that take as
% output an existing forecast run dataset and a provisional run dataset. It
% converts the structure datasets from MAPS to cell array datasets for 
% EASE.
%
% INPUTS:
%   -> Model: valid MAPS model
%   -> MAPSrunData: structure of forecast run data from MAPS
%   -> MAPSprovRunData: structure of provisional run data from MAPS
%
% OUTPUTS:
%   -> EASErunData: cell array of forecast run data for EASE
%   -> EASEprovRunData: cell array of provisional run data for EASE
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> convert_forecast_run_dataset_from_MAPS_to_EASE_format
%   -> generate_MAPS_exception_and_add_cause
%   -> convert_prov_run_dataset_from_MAPS_to_EASE_format
%
% DETAILS:
%   -> This helper translates the MAPS structure formats for forecast and
%      provisional run data to cell arrays for EASE.
%   -> See the contents of the two underlying functions for details.
%
% NOTES:
%   -> See xxxxxx for a discussion of EASE- & MAPS-represented forecast run
%      dataset formats in EASE and MAPS.
%   -> The underlying functions do not presume anything about which data
%      components must exist (i.e. they do not check that the content is
%      consistent with an ability to complete any particular operation, it
%      will just convert whatever it finds given the input and data
%      configuration information).
%
% This version: 21/02/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number & shape of inputs is as expected. All inpurs are
% compulsory. The EASE datasets must be cell arrays.
if nargin < 3
    errId = 'MAPS:convert_run_datasets_from_MAPS_to_EASE_format:BadNargin';
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = 'MAPS:convert_run_datasets_from_MAPS_to_EASE_format:BadInput1';
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(MAPSrunData)
    errId = 'MAPS:convert_run_datasets_from_MAPS_to_EASE_format:BadInput2';
    generate_and_throw_MAPS_exception(errId);
elseif ~isstruct(MAPSprovRunData)
    errId = 'MAPS:convert_run_datasets_from_MAPS_to_EASE_format:BadInput3';
    generate_and_throw_MAPS_exception(errId);
end

%% SETUP A MASTER EXCEPTION
% Setup a master exception to add all exception causes as encountered
% below.
masterErrId = ['MAPS:convert_run_datasets_from_MAPS_to_EASE_format',...
    ':BadRunData'];
BadRunDataE = generate_MAPS_exception(masterErrId);

%% CONVERT THE FORECAST RUN DATASET
% Attempt to convert the forecast run dataset input to EASE format. Capture
% any exceptions and add to the master exception.
try
    EASErunData = convert_forecast_run_dataset_from_MAPS_to_EASE_format(...
        Model,MAPSrunData);
catch DataE
    errId = [masterErrId,':Instance'];
    ForecastRunDataE = generate_MAPS_exception_and_add_cause(...
        DataE,errId,{'forecast'});
    BadRunDataE = addCause(BadRunDataE,ForecastRunDataE);
end

%% CONVERT THE PROVISIONAL RUN DATASET
% Attempt to convert the provsional run dataset input to EASE format.
% Capture any exceptions and add to the master exception.
try
    EASEprovRunData = convert_prov_run_dataset_from_MAPS_to_EASE_format(...
        Model,MAPSprovRunData);
catch DataE
    errId = [masterErrId,':Instance'];
    ProvRunDataE = generate_MAPS_exception_and_add_cause(...
        DataE,errId,{'Provisional'});
    BadRunDataE = addCause(BadRunDataE,ProvRunDataE);
end

%% THROW MASTER EXCEPTION
% Throw the master exception if any exceptions were encountered above.
if ~isempty(BadRunDataE.cause)
    throw(BadRunDataE);
end

end