function F = vou_to_pwcgc(A,V)

% Calculate time-domain pairwise-conditional Granger causality rates (Granger-
% causal graph) for a vector Ornstein-Uhlenbeck (VOU) process. Uses a state-space
% method which involves solving an associated continuous-time algebraic Riccati
% equation (CARE).
%
% A     - VOU coefficients matrix
% V     - VOU Wiener process covariance matrix
%
% F     - pairwise-conditional Granger causality rates (Granger-causal graph)
%
% REFERENCES:
%
% (1) L. Barnett and A. K. Seth (2015): Granger causality for state-space models, Phys. Rev. E 91(4) Rapid Communication.
% (2) L. Barnett and A. K. Seth (2016): Detectability of Granger causality for subsampled continuous-time neurophysiological processes, J. Neurosci. Methods 275.
% (3) L. Barnett (2017): Granger causality rate for a vector Ornstein-Uhlenbeck process (working notes).
%
% (C) Lionel Barnett, May 2017

[n, n1]  = size(A); assert(n1 == n, 'VOU coeffcicients matrix must be square');
[n1,n2]  = size(V); assert(n1 == n2,'VOU covariance matrix must be square');
                    assert(n1 == n, 'VOU covariance matrix must be same size as coefficients matrix');

F = nan(n);
for y = 1:n
    r = [1:y-1 y+1:n]; % omit y
	[P,~,~,rep] = icare(A(y,y)',A(r,y)',V(y,y),V(r,r),V(r,y)'); % NOTE: these are actually a quadratic equations for P!
	if vouerror(rep.Report), return; end % check CARE report, bail out on error
	F(r,y) = (A(r,y).^2).*(P./diag(V(r,r)));
end
