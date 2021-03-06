function [n,r,L] = ss_parms(A,C,K,V)

[r, r1] = size(A); assert(r1 == r,           'SS: bad ''A'' parameter');
[n, r1] = size(C); assert(r1 == r,           'SS: bad ''C'' parameter');
[r1,n1] = size(K); assert(n1 == n && r1 == r,'SS: bad ''K'' parameter');
if nargin > 3
	[n1,n2] = size(V); assert(n1 == n && n2 == n,'SS: bad ''V'' parameter');
	if nargout > 2
		[L,p] = chol(V,'lower');
		if p ~= 0, L = []; end % not positive-definite; caller to test
	end
end
