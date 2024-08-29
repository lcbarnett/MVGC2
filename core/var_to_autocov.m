% Autocovariance sequence for VAR model
%
% NOTES:
%
% Gamma_k is returned in G(:,:,k+1)
%
% Positive qmax is maximum lags; calculate sequence iteratively until
% within tolerance or maximum lags exceeded.
%
% Negative qmax is absolute number of lags.
%
% if qmax == 0, the covariance matrix Gamma_0 is returned in G.
%
% Default tolerance is machine fp precision relative to Gamma_0.

function [G,q] = var_to_autocov(A,V,qmax,tol)

if isvector(A)
	A = reshape(A,1,1,length(A));
end
[n,~,p] = size(A);

% Associated VAR(1) parameters

[A1,V1] = var_companion(A,V);

% Solve the Lyapunov equation for the covariance matrix of the associated VAR(1)

G1 = dlyap(A1,V1);

G = reshape(G1(1:n,:),n,n,p);

alags = qmax < 0; % -qmax is absolute number of lags
if alags, q = -qmax; else, q = qmax; end
q1 = q+1;

% We already have p-1 lags; if that's enough, truncate if necessary and return

if q < p
	G = G(:,:,1:q1);
	return
end

% Initialise reverse covariance sequence

R = zeros(p*n,n);
for k = 1:p
	R((k-1)*n+1:k*n,:) = G(:,:,p-k+1);
end

% Calculate autocovariances iteratively

A = reshape(A,n,p*n);
pn1 = (p-1)*n;

if alags % calculate recursively from p lags up to q lags

	G = cat(3,G,zeros(n,n,q1-p)); % pre-allocate
	for k = p+1:q1
		G(:,:,k) = A*R;
		R = [G(:,:,k);R(1:pn1,:)]; % update reverse covariance sequence
	end

else     % calculate recursively from p lags until convergence or maximum lags exceeded

	if nargin < 4 || isempty(tol), tol = eps; end
	k = p;
	isqrcov = 1./sqrt(diag(G(:,:,1)));                % inverse square roots of (zero-lag) covariances
	while max(abs(isqrcov.*G(:,:,k).*isqrcov')) > tol % maximum absolute autocorrelation
		if k > qmax
			fprintf(2,'WARNING: covariance sequence failed to converge (increase max. lags?)\n');
			q = k-1;
			return
		end
		k = k+1;
		G(:,:,k) = A*R;
		R = [G(:,:,k);R(1:pn1,:)]; % update reverse covariance sequence
	end
	q = k-1;

end
