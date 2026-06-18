function [X,E,mtrunc] = var_to_tsdata_res(A,res,m,N,mtrunc,decfac)

% A version of var_to_tsdata that takes a residuals time series res, rather than
% a residuals covariance matrix. Residuals are sampled uniformly randomly with
% replacement from res.

if isvector(A)
	A = A(:)';
	n = 1;
	p = length(A);
else
	[n,n1,p] = size(A);
	assert(n1 == n,'Bad VAR coefficients array');
end

[n1,k] = size(res);
assert(n1 == n && ismatrix(res) && k >= n,'Bad residuals matrix');

if nargin < 4 || isempty(N), N = 1; end % single trial

if nargin < 5 || isempty(mtrunc) % automatic calculation - transients decay with rate given by VAR spectral radius
    if nargin < 6 || isempty(decfac)
		decfac = 1;
	else
		assert(isscalar(decfac) && isnumeric(decfac) && decfac >= 0,'for automatic truncation estimation, decay factor must be a positive number');
	end
    mtrunc = decfac*var_decorrlen(A,res*res');
    assert(~isinf(mtrunc),'VAR unstable - can''t truncate!');
else
	assert(isscalar(mtrunc) && isint(mtrunc),'truncation parameter must be an integer');
	if mtrunc < 0 % calculate from autocovariance decay, with max lags = -mtrunc, and tolerance in decfac
		if nargin < 6
			decfac = [];
		else
			assert(isscalar(decfac) && isnumeric(decfac) && decfac >= 0,'for autocovariance truncation estimation, tolerance must be a positive number');
		end
		mtrunc = var_decorrlen(A,res*res',-mtrunc,[],decfac);
	end
end

if N > 1 % multi-trial
	E = zeros(n,mtrunc+m,N);
	for r = 1:N
		E(:,:,r) = res(:,randi(k,1,mtrunc+m)); % sample from supplied residuals with replacement
	end
	X = zeros(n,mtrunc+m,N);
	for r = 1:N
		X(:,:,r) = mvfilter([],A,E(:,:,r));
	end
else
	E = res(:,randi(k,1,mtrunc+m)); % sample from supplied residuals with replacement
	X = mvfilter([],A,E);
end

if mtrunc > 0
	X = X(:,mtrunc+1:mtrunc+m,:);
	if nargout > 1
		E = E(:,mtrunc+1:mtrunc+m,:);
	end
end
