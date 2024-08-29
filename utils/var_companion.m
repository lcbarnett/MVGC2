function [A1,V1] = var_companion(A,V)

if isvector(A)
	A = reshape(A,1,1,length(A));
end

[n,n1,p] = size(A);
assert(n1 == n,'Bad coefficients array');
pn1 = (p-1)*n;
A1 = [reshape(A,n,p*n); eye(pn1) zeros(pn1,n)]; % VAR coefficients for 1-lag problem

if nargout > 1
	assert(nargin > 1,'Must supply residuals covariance matrix');
	[n1,n2] = size(V);
	assert(n1 == n && n2 == n,'Bad residuals covariance matrix');
	V1 = [V zeros(n,pn1); zeros(pn1,n) zeros(pn1)];
end
