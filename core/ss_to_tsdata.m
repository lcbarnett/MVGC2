function [X,Z,E,mtrunc] = ss_to_tsdata(A,C,K,V,m,N,mtrunc,decfac)

if nargin < 6 || isempty(N), N = 1; end % single trial

if nargin < 7 || isempty(mtrunc) % automatic calculation - transients decay with rate given by VAR spectral radius
    if nargin < 8 || isempty(decfac)
		decfac = 1;
	else
		assert(isscalar(decfac) && isnumeric(decfac) && decfac >= 0,'for automatic truncation estimation, decay factor must be a positive number');
	end
    mtrunc = decfac*var_decorrlen(A,V);
    assert(~isinf(mtrunc),'unstable - can''t truncate!');
else
	assert(isscalar(mtrunc) && isint(mtrunc),'truncation parameter must be an integer');
	if mtrunc < 0 % calculate from autocovariance decay, with max lags = -mtrunc, and tolerance in decfac
		if nargin < 8
			decfac = [];
		else
			assert(isscalar(decfac) && isnumeric(decfac) && decfac >= 0,'for autocovariance truncation estimation, tolerance must be a positive number');
		end
		mtrunc = var_decorrlen(A,V,-mtrunc,[],decfac);
	end
end

[n,r,L] = ss_parms(A,C,K,V);

% Generate Gaussian innovations with covariance matrix V

mtot = m+mtrunc;
E = nan(n,mtot+1,N);
Z = nan(r,mtot  ,N);
X = nan(n,mtot  ,N);
for u = 1:N
    E(:,:,u) = L*randn(n,mtot+1);
    Z(:,:,u) = mvfilter([],A,K*E(:,1:mtot,u)); % state variable is AR(1) with innovations K*E (lagged)
    X(:,:,u) = C*Z(:,:,u) + E(:,2:mtot+1,u);   % observations variable uses unlagged innovations
end

E = E(:,2:mtot+1,:); % align innovations with X

if mtrunc > 0
    X = X(:,mtrunc+1:end,:);
    if nargout > 1
		Z = Z(:,mtrunc+1:end,:);
		if nargout > 2
			E = E(:,mtrunc+1:end,:);
		end
	end
end
