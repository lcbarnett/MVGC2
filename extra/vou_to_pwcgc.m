function F = vou_to_pwcgc(A,V)

% DESCRIPTION:
%
% Calculate time-domain pairwise-conditional Granger causality rates (Granger-
% causal graph) for a vector Ornstein-Uhlenbeck (VOU) process X(t):
%
%	dX(t) = AX(t)dt a dW(t)
%
% with autoregression matrix A, and W(t) a Wiener process with dW(t) ~ N(0,Vdt).
% NOTE: A need not be stable (i.e., may have eigenvalues with nonnegative real
% part).
%
% PARAMETERS:
%
% A     - VOU coefficients matrix
% V     - VOU Wiener process covariance matrix
%
% RETURN VALUE:
%
% F     - pairwise-conditional Granger causality rates (Granger-causal graph)
%
% REFERENCES:
%
% (1) L. Barnett and A. K. Seth (2015): Granger causality for state-space models, Phys. Rev. E 91(4) Rapid Communication.
% (2) L. Barnett and A. K. Seth (2016): Detectability of Granger causality for subsampled continuous-time neurophysiological processes, J. Neurosci. Methods 275.
% (3) L. Barnett (2017): Granger causality rate for a vector Ornstein-Uhlenbeck process (working notes).
%
% (C) Lionel Barnett, 2024

[n, n1]  = size(A); assert(n1 == n, 'VOU coefficients matrix must be square');
[n1,n2]  = size(V); assert(n1 == n2,'VOU covariance matrix must be square');
                    assert(n1 == n, 'VOU covariance matrix must be same size as coefficients matrix');

F = nan(n);
for y = 1:n
    r = [1:y-1 y+1:n]; % omit y

	if all(A(r,y) == 0)
		F(r,y) = 0;
		continue
	end

	L = chol(V(r,r));
	AOL = A(r,y)'/L;
	VOL = V(y,r)/L;
	a = AOL*AOL';
	b = AOL*VOL'-A(y,y);
	c = VOL*VOL'-V(y,y);
	P = (sqrt(b^2-a*c)-b)/a;

	F(r,y) = (A(r,y).^2).*(P./diag(V(r,r)));
end
