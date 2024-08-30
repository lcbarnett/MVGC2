function [m,rho] = var_decorrlen(A,V,maxlags,rho,tol)

% Calculate number of lags for autocorrelation to decay to negligible

if nargin < 3 || isempty(maxlags) % calculate using VAR spectral norm
	if nargin < 4 || isempty(rho), rho = specnorm(A); end
	if rho > 1-eps,
		m = Inf; % unstable!
	else
		if nargin < 5 || isempty(tol), tol = eps; end
		m = ceil((log(eps)-log(norm(V)))/log(rho));
		p = size(A,3);
		if m < p, m = p; end
	end
else                              % calculate using VAR autocovariance sequence (more reliable)
	assert(isscalar(maxlags) && isint(maxlags) && maxlags > 0,'maximum lags parameter must be a positive integer');
	if nargin < 5, tol = []; end
	[~,m] = var_to_autocov(A,V,maxlags,tol);
	if (nargin < 4 || isempty(rho)) && nargout > 1, rho = specnorm(A); end
end
