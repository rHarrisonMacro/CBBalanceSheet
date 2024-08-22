function y = grid2d(x)
% *************************************************************************
%
% GRID2D produces a two dimensional grid which references coordinates 
% across a multidimensional plane. Each column refers to a dimension on 
% that plane. A cell input of two vectors - a 25*1 and a 50*1
% - produces a matrix y of (25.50=100)*2. Note that there is also a C mex
% version of this coded up.
%
% INPUTS:  x - a 1*d cell array, where d is the number of vectors
% OUTPUTS: y - a prod(n)*d matrix, where prod(n) is the total number of
%              coordinates and d is the number of planes.
% 
% SEE ALSO: FUNDEFN, FUNNODES, FUNBASIS, FUNEVAL
%
% *************************************************************************

if ~iscell(x)
    error('input must be a cell array');
elseif size(x,1) < 1 || size(x,2) < 1
    error('cell must be non-empty');
elseif size(x,1) ~= 1
    error('cell must be arranged as a row vector');  
end
d = size(x,2);
n = zeros(1,d);
outputSigStr = '[';
for j = 1 : d
    xj = x{j}; 
    n(j) = size(xj,1);
    if size(xj,2) ~= 1
        error('cell elements must be arranged as col. vectors'); 
    end   
    outputSigStr = [outputSigStr,'grid',num2str(j)];                        %#ok<AGROW>
    if j < d
        outputSigStr = [outputSigStr,','];                                  %#ok<AGROW>
    else
        outputSigStr = [outputSigStr,']'];                                  %#ok<AGROW>
    end
end
y = zeros(prod(n),d);

% *************************************************************************

funcCallSig = [outputSigStr,'=ndgrid(x{:});'];
eval(funcCallSig);

% *************************************************************************

for j = 1 : d
    y(:,j) = eval(['grid',num2str(j),'(:);']);
end

% *************************************************************************

end