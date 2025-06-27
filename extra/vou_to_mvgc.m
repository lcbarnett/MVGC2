function [F,err] = vou_to_mvgc(A,V,x,y)

% DESCRIPTION:
%
% Calculate time-domain conditional Granger causality rate for a vector
% Ornstein-Uhlenbeck (VOU) process X(t):
%
%	dX(t) = AX(t)dt a dW(t)
%
% with autoregression matrix A, and W(t) a Wiener process with dW(t) ~ N(0,Vdt).
% NOTE: A need not be stable (i.e., may have eigenvalues with nonnegative real
% part).
%
% Source/target/conditioning variables may be multivariate. Uses a state-space
% method which involves solving an associated continuous-time algebraic Riccati
% equation (CARE).
%
% PARAMETERS:
%
% A     - VOU coefficients matrix
% V     - VOU Wiener process covariance matrix
% x     - multi-index of target variable
% y     - multi-index of source variable
%
% F     - Granger causality rate from y to x, conditional on other variables
% err   - CARE error report number (zero if no error)
%
% RETURN VALUES:
%
% Possible CARE errors are (see 'icare' in the Matlab Control System Toolbox):
%
%	0 - No errors (unique solution is accurate)
%	1 - Solution accuracy is poor
%	2 - Solution not finite
%	3 - No solution found (Hamiltonian spectrum has imaginery eigenvalues)
%	4 - No solution found ("pencil" is singular)
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

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),     'Some target indices out of range');
assert(all(y >=1 & y <= n),     'Some source indices out of range');
assert(isempty(intersect(x,y)), 'Source/target multi-indices overlap');

z = 1:n; z([x y]) = []; % indices of remaining (conditioning) variables
r = [x z];              % indices for the reduced system

err = 0;

if all(A(x,y) == 0)
	F = 0;
	return
end

if length(y) == 1 % P scalar, so CARE is a quadratic equation.

	L = chol(V(r,r));
	AOL = A(r,y)'/L;
	VOL = V(y,r)/L;
	a = AOL*AOL';
	b = AOL*VOL'-A(y,y);
	c = VOL*VOL'-V(y,y);
	P = (sqrt(b^2-a*c)-b)/a;
	F = P*trace(V(x,x)\(A(x,y)*A(x,y)'));

else

	[P,~,~,rep] = icare(A(y,y)',A(r,y)',V(y,y),V(r,r),V(r,y)');
	err = rep.Report;
	if err ~= 0
		F = NaN;
		return
	end
	F = trace(V(x,x)\(A(x,y)*P*A(x,y)'));

end
