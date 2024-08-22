function EASEdata = convert_dataset_from_MAPS_to_EASE_format(...
    Model,MAPSdata,MAPSdataConfig,EASEdataConfig)
% This helper converts data output from MAPS into EASE compatible format.
% It should be used as a helper in all EASE interfaces in MAPS that take
% provisional and forecast run dataset outputs from MAPS. It converts 
% the structure dataset from MAPS to a cell array dataset for EASE.
%
% INPUTS:
%   -> Model: valid MAPS model
%   -> MAPSdata: structure of data from MAPS
%   -> MAPSdataConfig: cell array of configuration info for the MAPS data 
%   -> EASEdataConfig: cell array of configuration info for the EASE data 
%
% OUTPUTS:
%   -> EASEdata: cell array of data for EASE
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_and_MAPS_exception
%   -> get_metadata_and_data_to_pack_in_EASE_data_cell (sub-function)
%   -> add_series_to_EASE_data_cell (sub-function)
%   -> generate_MAPS_exception_and_add_cause
%
% DETAILS:
%   -> This helper translates the standard MAPS structure format for 
%      datasets (including forecast run data and provisional run data) to a 
%      cell array for EASE. 
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
%   -> See xxxxxx for a discussion of EASE- & MAPS-represented forecast run 
%      dataset formats in EASE and MAPS.
%   -> This function does not presume anything about which data components 
%      must exist (i.e. it does not check that the content is consistent 
%      with an ability to complete any particular operation, it will just 
%      convert whatever it finds given the input and data configuration 
%      information). See the MAPS dataset content validation helper for a
%      function which does that.
%   -> It is straightforward to change details of the format of either the
%      MAPS- or EASE-formatted datasets by changing the configuration
%      information which is passed into this function as input (i.e. 
%      without touching this function). However, this function makes
%      several assumptions about the format of the configurations and the
%      general format of the data. For example, EASE datasets are always
%      assumed to be stored in cell arrays and MAPS in structures. Changing
%      that more general information would require major changes to this
%      function. See the code below for the assumptions made.
%
% This version: 14/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number & shape of inputs is as expected. All inpurs are 
% compulsory. The MAPS data must be a structure and the configurations must 
% be cell arrays with the MAPS data configuration containing at least 6 
% columns. 
if nargin < 4
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);      
elseif ~isstruct(MAPSdata)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif ~iscellstr(MAPSdataConfig) || ndims(MAPSdataConfig)~=2 || ...
        size(MAPSdataConfig,2)<6
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);     
elseif ~iscellstr(EASEdataConfig) || ndims(EASEdataConfig)~=2
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);     
end

%% COMPUTE THE NUMBER OF METADATA TAGS IN THE EASE DATA CONFIG
% The EASE data cell array must be configured so that the final column
% contains data, the penultimate column contains model metadata for that
% data and the first n columns contain metadata describing the data series.
% The configuration must include n columns of data series description to
% match the data.
nEASEdataIds = size(EASEdataConfig,2);

%% CHECK CONSISTENCY OF THE EASE AND MAPS DATA CONFIGURATION
% It must be possible to map each of the potential series in the MAPS data
% input into the EASE data output, so the EASE data configuration must have
% as many columns as the MAPS data configuration.
nPossDataTypes = size(MAPSdataConfig,1);
if size(EASEdataConfig,1) ~= nPossDataTypes
    errId =  ['MAPS:',mfilename,':EASEconfigInconsistentWithMAPSconfig'];
    generate_and_throw_MAPS_exception(errId);
end

%% CONSTRUCT EASE DATA
% Setup the EASE dataset cell and an exception to add any errors 
% encountered. Run through each of the possible MAPS data types. For each, 
% attempt to extract associated metadata and data. If any data is found, 
% call a sub-function to add the series to the MAPS structure. Catch any 
% errors encountered and add them to the master exception setup at the top. 
EASEdata = cell(0,5);
masterErrId = ['MAPS:',mfilename,':DataConversionFailure'];
DataConversionE = generate_MAPS_exception(masterErrId);
for iType = 1:nPossDataTypes
    try
        [iTypeMetadata,iTypeVecs] = ...
            get_metadata_and_data_to_pack_in_EASE_data_cell(...
            Model,MAPSdata,MAPSdataConfig(iType,:));
        if ~isempty(iTypeMetadata)
            EASEdata = add_series_to_EASE_data_cell(EASEdata,...
                iTypeMetadata,iTypeVecs,...
                EASEdataConfig(iType,1:nEASEdataIds));
        end
    catch iTypeConversionE
        errId = [masterErrId,':Instance'];
        iDataConversionE = generate_MAPS_exception_and_add_cause(...
            iTypeConversionE,errId,MAPSdataConfig(iType,1:3));
        DataConversionE = addCause(DataConversionE,iDataConversionE);
    end
end

%% THROW ANY EXCEPTION(S) ENCOUNTERED
% Throw the master exception created above if any data conversion 
% exceptions were encountered.
if ~isempty(DataConversionE.cause)
    throw(DataConversionE);
end

end

%% FUNCTION TO UNPACK DATA FROM MAPS STRUCTURE
function [seriesMetadata,seriesVecs] = ...
    get_metadata_and_data_to_pack_in_EASE_data_cell(...
    Model,MAPSdata,MAPSseriesConfig)
% This helper unpacks associated data series from the MAPS data structure.
% Together with the function below, it separates the code that adds the 
% data series to the EASE cell from the code that extracts the series from 
% the MAPS structure.
%
% INPUTS:
%   -> Model: valid MAPS model
%   -> MAPSdata: structure of data from MAPS
%   -> MAPSdataConfig: cell string array of config data for the MAPS data
%
% OUTPUTS:
%   -> seriesMetadata: cell string array of metadata for each series
%   -> seriesVecs: cell array of data vectors for each series
%
% CALLS:
%   -> unpack_data_from_MAPS_dataset
%   -> generate_and_throw_MAPS_exception
%   -> unpack_model

%% CALL THE MAPS DATASET UNPACKER HELPER
% Call the generic MAPS dataset structure unpacker to attempt to unpack the
% data described by the first three components of the MAPS dataset
% configuration. This function will return empty if the data does not exist
% within the dataset.
seriesData = unpack_data_from_MAPS_dataset(MAPSdata,MAPSseriesConfig(1:3));

%% SEPARATE THE DATA INTO NUMERIC AND METADATA
% Construct the data to add to the EASE cell. The method for construction
% will depend on the storage type string from the MAPS data configuration.
seriesMetadata = [];
seriesVecs = [];
if ~isempty(seriesData)
    seriesStorageType = MAPSseriesConfig{4};
    modelMetadataStr = MAPSseriesConfig{6};
    switch seriesStorageType
        case 'cell'
            if ~iscell(seriesData) || size(seriesData,2)~=2
                errId = ['MAPS:',mfilename,':BadCell'];
                generate_and_throw_MAPS_exception(errId);
            else
                seriesMetadata = seriesData(:,1);
                seriesVecs = seriesData(:,2);
            end
        case 'matrix'
            if ~isnumeric(seriesData) || ndims(seriesData)~=2
                errId = ['MAPS:',mfilename,':BadMatrix'];
                generate_and_throw_MAPS_exception(errId);
            else                
                seriesMetadata = unpack_model(Model,{modelMetadataStr});
                [nSeries,T] = size(seriesData);
                seriesVecs = mat2cell(seriesData,ones(1,nSeries),T);
            end
        otherwise
            errId = ['MAPS:',mfilename,':UnknownStorageType'];
            generate_and_throw_MAPS_exception(errId,{seriesStorageType});
    end
end

end

%% FUNCTION TO ADD SERIES TO THE EASE CELL DATASET
function EASEdata = add_series_to_EASE_data_cell(EASEdata,...
    seriesMetadata,seriesVecs,EASEseriesConfig)
% This helper adds associated data series to the EASE cell dataset.
% Together with the function above, it separates the code that adds the 
% data series to the EASE cell from the code that extracts the series from 
% the MAPS structure.
%
% INPUTS:
%   -> EASEdata: cell of data for EASE
%   -> seriesMetadata: cell string array of metadata for each series
%   -> seriesVecs: cell array of data vectors for each series
%   -> EASEseriesConfig: cell string array of config data for the EASE data
%
% OUTPUTS:
%   -> EASEdata: updated cell of data for EASE
%
% CALLS:
%   -> generate_and_throw_MAPS_exception

%% CHECK THE SERIES MNEMONICS AND VALUES
% Check that the series metadata and data have the correct data types.
if ~iscellstr(seriesMetadata)
    errId = ['MAPS:',mfilename,':BadMAPSmetadata'];
    generate_and_throw_MAPS_exception(errId);
elseif ~all(cellfun(@isnumeric,seriesVecs))
    errId = ['MAPS:',mfilename,':BadMAPSdata'];
    generate_and_throw_MAPS_exception(errId);
end

%% COMPUTE NUMBER OF SERIES TO ADD
% Determine the number of series to add. Check that the metadata that goes
% with the data has consistent dimensions.
nSeries = size(seriesVecs,1);
if size(seriesMetadata,1) ~= nSeries;
    errId = ['MAPS:',mfilename,':InconsistentMAPSmetadataData'];
    generate_and_throw_MAPS_exception(errId);
end

%% APPEND THE NEW DATA TO THE EASE DATA
% Append the new data to the EASE dataset. Expand the series metadata in
% the configuration input to cover the range of new data and add the
% metadata / data into the final two columns
EASEdata = [EASEdata;...
    [repmat(EASEseriesConfig,[nSeries 1]) seriesMetadata seriesVecs]];

end