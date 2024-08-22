function modelDataset = unpack_data_from_MAPS_dataset(MAPSdata,seriesIds)
% This helper unpacks data from a generic MAPS dataset structure.
% It is used as a helper in the conversion of MAPS datasets to EASE format
% and in the validation of the content of MAPS datasets in EASE interfaces.
%
% INPUTS:
%   -> MAPSdata: structure of data from MAPS
%   -> seriesIds: cell string array of series IDs for the data to unpack
%
% OUTPUTS:
%   -> modelDataset: dataset for a particular variable type
%
% CALLS:
%   -> generate_and_throw_MAPS_exception
%
% DETAILS:
%   -> This helper unpacks data associated with a particular model and
%      variable type in that model, as well as time horizon type (eg
%      forecast data for model variables).
%   -> MAPS represents the data in matrices (or cell arrays if the data is
%      not a complete model dataset) stored in a single structure.
%   -> The fieldnames is that structure describe what the data type is and
%      the series identifiers passed in to this function hierarchicly
%      describe fieldnames in the stucture and hence identify a particular
%      dataset.
%   -> For example, a complete dataset for the raw observables over the
%      past (ie back data) could be in the MAPS dataset in fields
%      'Past.rawObservables' in which case the first two elements of the
%      series identifiers passed in would be 'Past' & 'rawObservables' with
%      any subsequent elements being empty strings.
%   -> This function will return empty if the dataset does not exist in the
%      structure input.
%
% NOTES:
%   -> See xxxxxx for a discussion of EASE- & MAPS-represented forecast run
%      dataset formats.
%
% This version: 16/03/2011
% Author(s): Matt Waldron

%% CHECK INPUTS
% Check that the number & shape of inputs is as expected. All inpurs are 
% compulsory. The EASE data and configuration must be cell arrays and the 
% MAPS data configuration must be a four-column cell array. 
if nargin < 2
    errId = ['MAPS:',mfilename,':BadNargin'];
    generate_and_throw_MAPS_exception(errId,{num2str(nargin)});
elseif ~isstruct(MAPSdata)
    errId = ['MAPS:',mfilename,':BadInput1'];
    generate_and_throw_MAPS_exception(errId);      
elseif ~iscellstr(seriesIds) || size(seriesIds,1)~=1
    errId = ['MAPS:',mfilename,':BadInput2'];
    generate_and_throw_MAPS_exception(errId); 
end

%% SEARCH THROUGH MAPS DATASET STRUCTURE
% Search through the MAPS data structure for the model dataset. The search
% is made hierarchicly starting with the first layer of fields. If that
% field exists in the dataset input, then the data is unpacked from the 
% structure if the next series Id is empty (signifying that there are no 
% more sub-fields in that part of the structure) or, if not, the field is
% checked to see if it is a structure itself (because there are sub-fields
% in that part of the structure). At each point, if the field does not
% exist in the structure the model dataset output is set to empty.
MAPSdataToUnpack = MAPSdata;
nSeriesIds = size(seriesIds,2);
seriesIds = [seriesIds {''}];
for iId = 1:nSeriesIds
    if isfield(MAPSdataToUnpack,seriesIds{iId})
        if isempty(seriesIds{iId+1})
            modelDataset = MAPSdataToUnpack.(seriesIds{iId});
            break
        elseif ~isstruct(MAPSdataToUnpack.(seriesIds{iId}))
            errId = ['MAPS:',mfilename,':BadField'];
            generate_and_throw_MAPS_exception(errId,{seriesIds{iId}});
        else
            MAPSdataToUnpack = MAPSdataToUnpack.(seriesIds{iId});
        end
    else
        modelDataset = [];
        break
    end
end

end