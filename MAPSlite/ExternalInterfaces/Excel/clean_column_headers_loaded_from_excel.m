function colHeaders = clean_column_headers_loaded_from_excel(rawColHeaders)
% This helper cleans and checks column headers loaded in from Excel.
%
% INPUTS:
%   -> rawColHeaders: row cell array of data loaded as column header from
%      Excel
%
% OUTPUTS:
%   -> colHeaders: row cell array of cleaned and checked headers
%
% DETAILS:
%   -> This function is a helper used in MAPS' family of Excel loading
%      functions to clean and check column headers loaded in.
%   -> It cleans them by converting any NaN values to empty strings ('').
%      This is necessary because MATLAB reads empty cells in Excel in as
%      NaN and in an abitrary set of column headers, users may want to
%      allow some of those headers to be empty.
%   -> It then checks that the cleaned headers et constitutes a row cell
%      string array (which ensures that numerics are not part of the set of
%      headers loaded in).
%
% NOTES:
%   -> This is part of a family of functions in MAPS which load in and save
%      out data.
%
% This version: 16/01/2013
% Author(s): Matt Waldron

%% CHECK INPUTS
if nargin < 1
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~is_row_cell_array(rawColHeaders)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);
end

%% SEARCH FOR EMPTY ELEMENTS REPLACE WITH STRING EQUIVALENT
% Empty cells in Excel are read in a NaN by MATLAB, so we are looking for
% numeric NaNs and these are replaced by ''.
isNumericColHeader = cellfun(@(x) isnumeric(x),rawColHeaders);
isEmptyColHeader = isNumericColHeader;
isEmptyColHeader(isNumericColHeader) = ...
    cellfun(@(x) isnan(x),rawColHeaders(isNumericColHeader));
colHeaders = rawColHeaders;
colHeaders(isEmptyColHeader) = {''};

%% CHECK COLUMN HEADERS ARE VALID (STRINGS)
if ~is_row_cell_string_array(colHeaders)
    errId = ['MAPS:',mfilename,':BadColHeaders'];
    generate_and_throw_MAPS_exception(errId);
end

end