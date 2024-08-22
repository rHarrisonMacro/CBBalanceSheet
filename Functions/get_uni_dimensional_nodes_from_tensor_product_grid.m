function snodes = get_uni_dimensional_nodes_from_tensor_product_grid(s)
% Extracts uni-dimensional grids from tensor product formed multi grid.
%
% INPUTS:   
%   -> s: ns*ds multi-dimensional grid matrix formed by tensor product
%
% OUTPUTS:  
%   -> snodes: 1*ds cell array of uni-dimensional nodes

%% INITIALISE OUTPUT
ds = size(s,2);
snodes = cell(1,ds);

%% GET UNIQUE NODES BY DIMENSION
for is = 1:ds
    snodes{is} = unique(s(:,is));
end

end