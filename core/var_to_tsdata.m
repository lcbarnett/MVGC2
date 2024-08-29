%% var_to_tsdata
%
% Generate random multi-trial Gaussian VAR time series
%
% <matlab:open('var_to_tsdata.m') code>
%
%% Syntax
%
%     [X,E,mtrunc] = var_to_tsdata(A,V,m,N,mtrunc,decfac)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     V          residuals covariance matrix
%     m          number of observations per trial
%     N          number of trials (default: 1)
%     mtrunc     number of initial time observations to truncate or (default) empty for automatic calculation
%     decfac     initial transients decay factor (default: 1)
%
% _output_
%
%     X          multi-trial Gaussian VAR time series
%     E          residuals time series
%     mtrunc     actual number of initial time steps truncated
%
%% Description
%
% Return |N| time series of length |m| sampled from a VAR model with
% coefficients matrix |A|, and iid Gaussian residuals with covariance matrix
% |V|:
%
% <<eq_var.png>>
%
% If mtrunc is supplied it is taken to be the the number of initial
% (non-stationary transient) observations to truncate; otherwise (default) the
% spectral radius of A (see function <specnorm.html |specnorm|>) is
% calculated and used to estimate a suitable number mtrunc of observations to
% assumed stationarity (roughly till autocorrelation decays to its stationary
% value); set decfac > 1 for longer equilibriation. If mtrunc is negative, then
% autocovariance decay is utilised, with maximum lags = -mtrunc, and tolerance
% in decfac. If mtrunc takes the special value 'stationary', the first p values
% of the time series are generated from the companion VAR(1), and then truncated.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
%% See also
%
% <specnorm.html |specnorm|> |
% <mvfilter.html |mvfilter|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [X,E,mtrunc] = var_to_tsdata(A,V,m,N,mtrunc,decfac)

if isvector(A)
	A = A(:)';
	n = 1;
	p = length(A);
else
	[n,n1,p] = size(A);
	assert(n1 == n,'Bad VAR coefficients array');
end

if nargin < 4 || isempty(N), N = 1; end % single trial

statts = false;
if nargin < 5 || isempty(mtrunc) % automatic calculation - transients decay with rate given by VAR spectral radius
    if nargin < 6 || isempty(decfac)
		decfac = 1;
	else
		assert(isscalar(decfac) && isnumeric(decfac) && decfac >= 0,'for automatic truncation estimation, decay factor must be a positive number');
	end
    mtrunc = decfac*var_decorrlen(A,V);
    assert(~isinf(mtrunc),'VAR unstable - can''t truncate!');
else
	if ischar(mtrunc)
		assert(strcmpi(mtrunc,'stationary'),'unknown truncation parameter');
		statts = true;
		mtrunc = p;
	else
		assert(isscalar(mtrunc) && isint(mtrunc),'truncation parameter must be an integer');
		if mtrunc < 0 % calculate from autocovariance decay, with max lags = -mtrunc, and tolerance in decfac
			if nargin < 6
				decfac = [];
			else
				assert(isscalar(decfac) && isnumeric(decfac) && decfac >= 0,'for autocovariance truncation estimation, tolerance must be a positive number');
			end
			mtrunc = var_decorrlen(A,V,-mtrunc,[],decfac);
		end
	end
end

[VL,cholp] = chol(V,'lower');
assert(cholp == 0,'covariance matrix not positive-definite');

if statts % generate first p stationary observations
	[A1,V1] = var_companion(A,V); % associated VAR(1) parameters
	G1 = dlyap(A1,V1);            % Solve the Lyapunov equation for the covariance matrix of the associated VAR(1)
	[G1L,cholp] = chol(G1,'lower');
	assert(cholp == 0,'VAR(1) covariance matrix not positive-definite');
	pn = p*n;
	if N > 1
		E = zeros(n,p+m,N);
		for r = 1:N
			E(:,:,r) = [reshape(G1L*randn(pn,1),n,p) VL*randn(n,m)];
		end
	else
		E = [reshape(G1L*randn(pn,1),n,p) VL*randn(n,m)];
	end
else
	if N > 1 % multi-trial
		E = zeros(n,mtrunc+m,N);
		for r = 1:N
			E(:,:,r) = VL*randn(n,mtrunc+m);
		end
	else
		E = VL*randn(n,mtrunc+m);
	end
end

if N > 1 % multi-trial
	X = zeros(n,mtrunc+m,N);
	for r = 1:N
		X(:,:,r) = mvfilter([],A,E(:,:,r));
	end
else
	X = mvfilter([],A,E);
end

if mtrunc > 0
	X = X(:,mtrunc+1:mtrunc+m,:);
	if nargout > 1
		E = E(:,mtrunc+1:mtrunc+m,:);
	end
end
