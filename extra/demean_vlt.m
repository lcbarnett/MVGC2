function [X,xmean,xstd] = demean_vlt(X,normalise)

% Temporal de-mean and optionally normalise by variance for variable-length trials data

if nargin < 2 || isempty(normalise), normalise = false; end

assert(isvector(X) && iscell(X),'X must be a cell vector of matrices');
n = size(X{1},1);
N = length(X);
m = zeros(N,1);
for r = 1:N
	assert(ismatrix(X{r}),'trials data must be matrices');
	[n1,m(r)] = size(X{r});
	assert(n1 == n,'number of variables in trials don''t match');
end

o = [0;cumsum(m)]; % trial length offsets
M = o(end);

Y = zeros(n,M);
for r = 1:N
	Y(:,o(r)+1:o(r+1)) = X{r};
end
xmean = mean(Y,2);
Y  = bsxfun(@minus,Y,xmean);

if normalise
	xstd = std(Y,[],2);
	Y = bsxfun(@rdivide,Y,xstd);
elseif nargout > 2
	xstd = std(Y,[],2);
end

for r = 1:N
	X{r} = Y(:,o(r)+1:o(r+1));
end
