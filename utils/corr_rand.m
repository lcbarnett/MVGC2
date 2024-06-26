function [R,L,grerr,retries,iters] = corr_rand(n,g,vexp,tol,maxretries)

% Generate random correlation matrix with specified multi-information
% ("log-generalised correlation" - see multiinfo.m). Creates a uniform random
% orthogonal matrix and a random variance vector with independent chi^2(1)
% distribution, then adjusts the variance vector until within tolerance of
% specified multi-information.
%
% The multi-information of a correlation matrix R is g = -log|R|. For n = 2,
% g = -log(1-rho^2) where rho is the Pearson correlation coefficient (so for
% small rho, g ~ rho^2.
%
% If g is not supplied (or is empty) a correlation matrix is sampled uniformly
% at random over the space of n x n correlation matrices using the "onion"
% algorithm. If g is zero, the identity matrix is returned. If g is negative,
% then -g is taken as a factor applied to the mean multi-information of a
% uniform random sampling of the space of correlation matrices.
%
% NOTE: if you want a covariance matrix with standard deviations in the (column)
% vector s, then use V = s.*R.*s', or L -> s.*L
%
% n          - number of dimensions
% g          - multi-information (g = -log|R|); g = 0 implies zero correlation (default: empty)
% vexp       - variance exponent (see 'corr_rand_exponent'; default: 2)
% tol        - numerical tolerance (default: sqrt(eps))
% maxretries - maximum retries to find large enough correlation (default: 1000)
%
% R          - correlation matrix
% L          - Cholesky (left) factor: L*L' = R
% grerr      - relative error of actual multi-information
% retries    - number of retries
% iters      - number of binary chop iterations

if nargin < 2                         g          = [];        end
if nargin < 3 || isempty(vexp),       vexp       = 2;         end
if nargin < 4 || isempty(tol),        tol        = sqrt(eps); end
if nargin < 5 || isempty(maxretries), maxretries = 1000;      end

retries = [];
iters   = [];
grerr   = [];

if isempty(g)
	R = onion(n); % uniform random on manifold of n x n correlation matrices (see 'onion.m')
	[L,pchol] = chol(R,'lower');
	if pchol ~= 0
		fprintf(2,'ERROR: ''corr_rand'' result not positive-definite (onion failed)\n');
		R = NaN;
		L = NaN;
		return
	end
	R = L*L';
	return
end

assert(isscalar(g) && isnumeric(g),'(log-generalised) correlation must be empty or a non-negative scalar');

if g < eps && g > -eps % effectively zero
	R = eye(n);
	L = eye(n);
	return
end

if g < 0 % g is a factor applied to the mean multi-information of a uniform sampling of correlation matrices
	g = -g*onion(n,[],true);
end

% We calcluate a (positive-definite) covariance matrix with given generalised
% correlation, and finally convert it to a correlation matrix. If V is a
% variance-covariance matrix, then
%
% g = sum(log(diag(V)))-logdet(V))

gtarget = g; % target value for g

% Find random orthogonal M and vector v of variances such that for V = M*diag(v)*M'
% (which will be pos-def), g is >= gtarget. The rationale is that adding the same
% constant c to all variances always decreases g, so we may use a binary chop to
% home in on gtarget. Note that since M is orthogonal, |V| = prod(v), so that we
% may calculate g efficiently as sum(log(diag(V))-log(v)).

for retries = 0:maxretries
	[Q,R] = qr(randn(n));
	v = realpow(abs(randn(n,1)),vexp);
	M = Q*diag(sign(diag(R))); % M orthogonal
	V = M*diag(v)*M';          % V is pos-def
	g = sum(log(diag(V))-log(v));
	if g >= gtarget
		break
	end
end
if g < gtarget
    fprintf(2,'ERROR: ''corr_rand'' timed out on retries (g too large?)\n');
    L = NaN;
    R = NaN;
    return
end
D = diag(V);

% Got M and v (for the binary chop we just need the variances D of V and the
% unrotated variances v). Now set binary chop initial high value so that g <
% gtarget; start c at 1 and keep doubling

iters = 0;
c = 1;
while g > gtarget
	g = sum(log(D+c) - log(v+c));
	iters = iters+1;
	c = 2*c;
end

% Do binary chop

chi = c;
clo = 0;
while true % binary chop
	c = (clo+chi)/2;
	g = sum(log(D+c) - log(v+c));
	iters = iters+1;
	if     g < gtarget % g too small, set chi to current
		chi = c;
	elseif g > gtarget % g too big, set clo to current
		clo = c;
	else
		break;
	end
	if chi-clo < tol
		break % run out of wiggle-room!
	end
end
% We have a c that is good as it gets

grerr = abs(gtarget-g)/gtarget; % relative error for multi-information

V = M*diag(v+c)*M';

% Check V really is pos-def

[L,pchol] = chol(V,'lower');
if pchol ~= 0
	fprintf(2,'ERROR: ''corr_rand'' result not positive-definite\n');
	R = NaN;
	L = NaN;
	return
end

% Convert to correlation matrix

L = diag(1./sqrt(sum(L.*L,2)))*L;
R = L*L';
