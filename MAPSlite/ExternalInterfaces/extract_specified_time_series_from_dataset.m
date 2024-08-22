function [dataMat,dataHeaders] = ...
    extract_specified_time_series_from_dataset(...
    allDataMat,allDataHeaders,seriesIds)
% This helper extracts specified time series from a set of time series.
%
% INPUTS:
%   -> allDataMat: nVars*T matrix of time series data
%   -> allDataHeaders: nVars*1 cell arry of identifiers for those series
%   -> seriesIds: string or cell string array of identifiers for the series
%      to extract
%
% OUTPUTS:
%   -> dataMat: nVars*T matrix of time series data
%   -> dataHeaders: nVars*1 column cell string array of headers consistent
%      with the seriesIds
%
% DETAILS:
%   -> This function extracts a specified subset of time series data from 
%      a MAPS time series dataset.
%   -> It searches for the specified series among the full set of series
%      and then extracts them in the order requested.
%
% NOTES:
%   -> This is part of a family of functions in MAPS which manage time 
%      series data.
%
% This version: 05/11/2012
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 3
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~is_real_two_dimensional_numeric_matrix(allDataMat)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
elseif ~is_column_cell_string_array(allDataHeaders) || ...
        ~isempty(find_repeated_strings(allDataHeaders))
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId);   
elseif size(allDataMat,1) ~= size(allDataHeaders,1)
    errId = ['MAPS:',mfilename,':InconsistentDataHeadersDims'];
    generate_and_throw_MAPS_exception(errId);     
elseif ~ischar(seriesIds) && ...
        ~is_column_cell_string_array(seriesIds) && ...
        ~is_row_cell_string_array(seriesIds)
    errId = ['MAPS:',mfilename,':BadInput3'];
    generate_and_throw_MAPS_exception(errId);   
end

%% DEAL WITH ALTERNATIVE SHAPE OF SERIES IDS
if ischar(seriesIds)
    dataHeaders = {seriesIds};
elseif size(seriesIds,1) == 1
    dataHeaders = seriesIds';
else
    dataHeaders = seriesIds;
end

%% ATTEMPT TO FIND INDICES OF THE DATA REQUESTED
try
    dataInds = lookup_index_numbers_in_string_array(...
        allDataHeaders,dataHeaders);
catch LookupE
    errId = ['MAPS:',mfilename,':UnknownIdentifiers'];    
    generate_MAPS_exception_add_cause_and_throw(LookupE,errId); 
end

%% EXTRACT RELEVANT DATA SERIES IN THE CORRECT ORDER
dataMat = allDataMat(dataInds,:);

end