function [A,V,E] = tsdata_to_var_vlt(X,p)

% Estimate VAR model for variable-length trials data (OLS only)

assert(isvector(X) && iscell(X),'X must be a cell vector of matrices');
n = size(X{1},1);
N = length(X);
m = zeros(N,1);
for r = 1:N
	assert(ismatrix(X{r}),'trials data must be matrices');
	[n1,m(r)] = size(X{r});
	assert(n1 == n,'number of variables in trials don''t match');
end
assert(all(p < m),'too many lags');

X = demean_mt(X);              % remove temporal mean

% store lags

p1 = p+1;
o = [0;cumsum(m-p)];           % trial length - model order offsets
M = o(end);
X0 = zeros(n,M);
for r = 1:N                    % concatenate trials for unlagged observations
	X0(:,o(r)+1:o(r+1)) = X{r}(:,p1:m(r));
end
XL = zeros(n,p,M);
for k = 1:p                    % concatenate trials for k-lagged observations
	for r = 1:N
		XL(:,k,o(r)+1:o(r+1)) = X{r}(:,p1-k:m(r)-k);
	end
end
XL = reshape(XL,n*p,M);        % stack lags

% OLS

A = X0/XL;                     % OLS (via QR decomposition)
if ~all(isfinite(A(:)))        % something went badly wrong
	V = [];
	E = [];
	return
end

if nargout > 1                 % calculate residuals covariance matrix
	EE  = X0-A*XL;             % residuals
	V = (EE*EE')/(M-1);        % residuals covariance matrix
	if nargout > 2             % calculate residuals
		E = cell(size(X));     % put residuals back into per-trial form
		for r = 1:N
			E{r} = EE(:,o(r)+1:o(r+1));
		end
	end
end

A = reshape(A,n,n,p);          % so A(:,:,k) is the k-lag coefficients matrix
