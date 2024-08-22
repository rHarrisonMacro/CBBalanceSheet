function MAPSdata = convert_dataset_from_EASE_to_MAPS_format(...
    Model,EASEdata,EASEdataConfig,MAPSdataConfig)
% This helper converts data input from EASE into a MAPS compatible format.
% It converts a standard cell array dataset from EASE to a standard 
% structure dataset for MAPS.
%
% INPUTS:
%   -> Model: valid MAPS model
%   -> EASEdata: cell array of data from EASE
%   -> EASEdataConfig: cell array of configuration info for the EASE data 
%   -> MAPSdataConfig: cell array of configuration info for the MAPS data 
%
% OUTPUTS:
%   -> MAPSdata: structure of data for MAPS
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%   -> generate_MAPS_exception
%   -> cell_array_vlookup
%   -> add_series_to_MAPS_data_structure (sub-function)
%   -> generate_MAPS_exception_and_add_cause
%
% DETAILS:
%   -> This helper translates the standard EASE cell array format for 
%      datasets (including forecast and provisional run datasets) to a 
%      structure for use in MAPS.
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
%      dataset formats.
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
% compulsory. The EASE data and configuration must be cell arrays and the 
% MAPS data configuration must be a four-column cell array. 
if nargin < 4
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(Model)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);      
elseif ~iscell(EASEdata)
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);
elseif ~iscellstr(EASEdataConfig) || ndims(EASEdataConfig)~=2
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);
elseif ~iscellstr(MAPSdataConfig) || ndims(MAPSdataConfig)~=2 || ...
        size(MAPSdataConfig,2)<6
    errId = ['MAPS:',mfilename,':BadInput4'];
    generate_and_throw_MAPS_exception(errId);      
end

%% CHECK CONSISTENCY OF EASA DATA WITH THE CONFIGURATION
% The EASE data cell array must be configured so that the final column
% contains data, the penultimate column contains model metadata for that
% data and the first n columns contain metadata describing the data series.
% The configuration must include n columns of data series description to
% match the data.
nEASEdataIds = size(EASEdataConfig,2);
if size(EASEdata,2) ~= nEASEdataIds+2
    errId =  ['MAPS:',mfilename,':EASEdataInconsistentWithConfig'];
    generate_and_throw_MAPS_exception(errId);
end

%% CHECK CONSISTENCY OF THE EASE AND MAPS DATA CONFIGURATION
% It must be possible to map each of the series in the EASE data input into
% the MAPS data output, so the MAPS data configuration must have as many 
% columns as the EASE data configuration.
nPossDataTypes = size(EASEdataConfig,1);
if size(MAPSdataConfig) ~= nPossDataTypes
    errId =  ['MAPS:',mfilename,':EASEconfigInconsistentWithMAPSconfig'];
    generate_and_throw_MAPS_exception(errId);
end

%% CONSTRUCT MAPS DATA
% Setup the MAPS structure and an exception to add any errors encountered. 
% Run through each of the possible EASE data types. For each, attempt to 
% extract associated metadata and data. If any data is found, call a 
% sub-function to add the series to the MAPS structure. Catch any errors 
% encountered and add them to the master exception setup at the top. 
MAPSdata = struct;
masterErrId = ['MAPS:',mfilename,':DataConversionFailure'];
DataConversionE = generate_MAPS_exception(masterErrId);
for iType = 1:nPossDataTypes
    try
        [iTypeMetadata,iTypeVecs] = cell_array_vlookup(...
            EASEdataConfig(iType,1:nEASEdataIds),...
            EASEdata(:,1:nEASEdataIds),...
            EASEdata(:,nEASEdataIds+1),EASEdata(:,nEASEdataIds+2));
        if ~isempty(iTypeMetadata)
            MAPSdata = add_series_to_MAPS_data_structure(...
                MAPSdata,Model,iTypeMetadata,iTypeVecs,...
                MAPSdataConfig(iType,:));
        end
    catch iTypeConversionE
        errId = [masterErrId,':Instance'];
        iDataConversionE = generate_MAPS_exception_and_add_cause(...
            iTypeConversionE,errId,EASEdataConfig(iType,1:nEASEdataIds));
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

%% FUNCTION TO ADD AN ASSOCIATED SET OF DATA SERIES TO MAPS STRUCTURE
function MAPSdata = add_series_to_MAPS_data_structure(MAPSdata,...
    Model,seriesMetadata,seriesVecs,MAPSdataConfig)
% This helper adds associated data series to the MAPS data structure.
% It separates the code that adds the data series to the MAPS structure
% from the code that extracts the series from the EASE data cell array.
%
% INPUTS:
%   -> MAPSdata: structure of data for MAPS
%   -> Model: valid MAPS model
%   -> seriesMetadata: cell string array of metadata for each series
%   -> seriesVecs: cell array of data vectors for each series
%   -> MAPSdataConfig: cell string array of MAPS config information
%
% OUTPUTS:
%   -> MAPSdata: updated structure of data for MAPS
%
% CALLS:
%   -> unpack_model
%   -> generate_and_throw_MAPS_exception
%   -> convert_data_cell_vectors_to_model_ordered_matrix (sub-function)

%% CHECK THE SERIES MNEMONICS AND VALUES
% Check that the series metadata and data have the correct data types.
if ~iscellstr(seriesMetadata)
    errId = ['MAPS:',mfilename,':BadEASEmetadata'];
    generate_and_throw_MAPS_exception(errId);
elseif ~all(cellfun(@isnumeric,seriesVecs))
    errId = ['MAPS:',mfilename,':BadEASEdata'];
    generate_and_throw_MAPS_exception(errId);  
end

%% UNPACK THE MAPS DATA CONFIGURATION
% Unpack all of the components of the MAPS configuration for the series
% passed in.
seriesType = MAPSdataConfig{1};
seriesSubType = MAPSdataConfig{2};
seriesSubSubType = MAPSdataConfig{3};
seriesStorageType = MAPSdataConfig{4};
modelMetadataStr = MAPSdataConfig{6};

%% CONSTRUCT THE DATA TO BE STORED
% Construct the data to add to the MAPS structure, whose format depends on
% the storage type string from the MAPS data configuration.
switch seriesStorageType
    case 'cell'
        seriesData = [seriesMetadata seriesVecs];
    case 'matrix'
        modelMetadata = unpack_model(Model,{modelMetadataStr});
        seriesData = convert_data_cell_vectors_to_model_ordered_matrix(...
            seriesMetadata,seriesVecs,modelMetadata);
    otherwise
        errId = ['MAPS:',mfilename,':UnknownStorageType'];
        generate_and_throw_MAPS_exception(errId,{seriesStorageType});
end

%% ADD THE DATA TO THE MAPS STRUCTURE
% Add the series to the MAPS data structure, taking into account variations
% in the potential depth of the structure.
if ~isempty(seriesSubSubType)
    MAPSdata.(seriesType).(seriesSubType).(seriesSubSubType) = seriesData;
elseif ~isempty(seriesSubType)
    MAPSdata.(seriesType).(seriesSubType) = seriesData;
else
    MAPSdata.(seriesType) = seriesData;
end

end
  
%% FUNCTION TO CONVERT DATA CELL VECTORS TO A MODEL ORDERED MATRIX
function dataMatrix = convert_data_cell_vectors_to_model_ordered_matrix(...
    seriesMetadata,seriesVecs,modelMetadata)
% This helper converts a cell array of data vectors to a data matrix.
% It reorders the vectors to put them in the same order as the series
% appear in the model. It then converts the cell array of vectors to a
% matrix.
%
% INPUTS:
%   -> seriesMetadata: cell array of metadata for each data series
%   -> seriesVecs: cell array of data series
%   -> modelMetadata: cell string array of metadata for each series
%
% OUTPUTS:
%   -> dataMatrix: matrix of data series
%
% CALLS:
%   -> lookup_model_index_numbers
%   -> generate_MAPS_exception_add_cause_and_throw

%% REORDER THE VECTORS OF VALUES TO BE CONSISTENT WITH THE MODEL
% Call a helper to find the index numbers of each of the variables 
% extracted from the data cell array from EASE in the model being used. 
% Reorder the vectors of values to match the ordering in the model. Note 
% that the aim is to reorder the unsorted cell array of values to match the 
% order in the model (rather than to extract information from the model) so 
% the inputs to the index helper are the reverse of what they are in most 
% index lookup calls.
try
    dataInds = lookup_model_index_numbers(seriesMetadata,modelMetadata);
    dataVecsSorted = seriesVecs(dataInds);
catch IndexLookupE
    errId = ['MAPS:',mfilename,':MetadataModelMismatch'];
    generate_MAPS_exception_add_cause_and_throw(IndexLookupE,errId);
end

%% CONVERT CELL ARRAYS TO MATRIX
% Convert the sorted cell array of vectors into a matrix using the MATLAB
% cell2mat command.
try
    dataMatrix = cell2mat(dataVecsSorted);
catch Cell2matE
    errId = ['MAPS:',mfilename,':CellToMatConversionFailure'];
    generate_MAPS_exception_add_cause_and_throw(Cell2matE,errId);
end  

end