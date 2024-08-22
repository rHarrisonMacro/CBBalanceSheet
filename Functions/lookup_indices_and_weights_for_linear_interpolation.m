function [inds,weights] = ...
    lookup_indices_and_weights_for_linear_interpolation(sLookup,snodes)
% Looks-up index numbers and weights for multi-dimensional linear interp.
% Also valid for linear extrapolation.  In multi-dimensional form, it is
% only valid under tensor product schemes.
%
% INPUTS:   
%   -> sLookup: 1*ds vector of function arguments to look up
%   -> snodes: 1*ds cell array of uni-dimensional nodes
%
% OUTPUTS:  
%   -> inds: ds^2 vector of index numbers required for the interp/extrap
%   -> weights: ds^2 vector of weights to apply to values at those indices

%% INITIALISE COLUMN-BY-COLUMN CELLS
ds = size(sLookup,2);
indsByCol = cell(1,ds);
weightsByCol = cell(1,ds);

%% COMPUTE INDICES AND WEIGHTS IN EACH COLUMN
nscum = 1;
for is = 1:ds
    isnodes = snodes{is};
    isinds = lookup_indices_from_table(isnodes,sLookup(is),3);
    isweights = (isnodes(isinds+1)-sLookup(is))/...
        (isnodes(isinds+1)-isnodes(isinds));
    indsByCol{is} = nscum*[isinds-1;isinds];
    weightsByCol{is} = [isweights;1-isweights];
    nscum = nscum*size(isnodes,1);
end

%% COMBINE INDICES AND WEIGHTS TOGETHER
inds = sum(grid2d(indsByCol),2)+1;
weights = prod(grid2d(weightsByCol),2);

end