function ind = lookup_indices_from_table(tabvals,x,endadj)
% LOOKUP performs a table lookup and outputs the following:
%           ind(i,j) = max(k) such that x(i,j) >= tabvals(k)
% 
% INPUTS:  tabvals - a sorted table vector of n values
%          x - an array of values
%          endadj - endpoint adjustment
% OUTPUTS: ind - an array of indices (the same size as x)
%
% ind = lookup(tabvals,x) does both endpoint adjustments as default 
%
% ind = lookup(tabvals,x,endadj) does no end point adjustment if endadj is 
% 0, adjusts the bottom end so that values of x < min(tabvals) are given an 
% index of 1 if endadj is 1, adjusts the top end so that values of 
% x > max(tabvals) are given an index of n-length(tabvals==tabvals(end)) if
% endadj is 2, and does both endpoint adjustment 1 and 2 if endadj = 3.
%
% The output of LOOKUP can be used to find the nearest table value as
% follows:
%           ind = ind+(x-tabvals(ind)>tabvals(ind+1)-x);
%           nearest = tabvals(ind);

%% CHECK INPUTS
if nargin < 2 || nargin > 3
    error('incorrect number of inputs');
elseif ~isnumeric(tabvals) || ~isnumeric(x)     
    error('tabvals and x need to be numeric matrices');
elseif ~isvector(tabvals) 
    error('tabvals must be a vector');
elseif any(diff(tabvals)<0)
    error('tabvals must be sorted in ascending order');
end

%% DEFINE DEFAULT FOR END ADJUSTMENT
if nargin == 2; 
    endadj = 3; 
end

%% COMPUTE INDICES
n = length(tabvals);
m = numel(x);
[temp,ind] = sort([tabvals(:);x(:)]);
temp = find(ind>n);
j = ind(temp)-n;
ind = reshape(temp-(1:m)',size(x));
ind(j) = ind(:);

%% ADJUST END POINTS
if endadj == 1 || endadj == 3
    ind(ind==0) = length(find(tabvals==tabvals(1))); 
end
if endadj == 2 || endadj == 3
    ind(ind==n) = n-length(find(tabvals==tabvals(n)));
end