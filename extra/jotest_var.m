function stats = jotest_var(A,V,normevs)

% Calculate eigenvalues for VECM likelihood function from VAR parameters.

if nargin < 3 || isempty(normevs), normevs = true; end % normalise eigenvectors?

[n,~,p] = size(A);

G = var_to_autocov(A,V,-p);    % n x n x (p+1) : G_0,..., G_p (-p means return exactly p lags)

P = G(:,:,1:p) - G(:,:,2:p+1); % n x n x p     : P_0,..., P_{p-1}

Q = zeros(n,n,p);              % n x n x p     : Q_0,..., Q_{p-1}
Q(:,:,1)   = P(:,:,1) + P(:,:,1)'; % Q_0 = P_0 + P_0'
Q(:,:,2:p) = P(:,:,2:p) - P(:,:,1:p-1);

% Set up regressions
%
% NOTE: For a potential method which does not have to store a full block-Toeplitz
% covariance matrix, see experimental/btsolve.m in https://github.com/lcbarnett/MVGC2
% which uses an algorithm described in Hirotugu Akaike, "Block Toeplitz Matrix Inversion",
% SIAM J. Appl. Math. 24(2):234-241, 1973.

QBT = sbtoeplitz(Q(:,:,1:p-1)); % symmetric block-Toeplitz of Q_0, ..., Q_{p-2}

QBTCHOL = chol(QBT); % upper Cholesky factor, so QBTCHOL'*QBTCHOL = QBT

np1 = (p-1)*n;
J = reshape(Q(:,:,2:p  ),n,np1)/QBTCHOL;
K = reshape(P(:,:,1:p-1),n,np1)/QBTCHOL;

% (Cross-)covariance matrices

S00 =  Q(:,:,1) - J*J';
S01 = -P(:,:,1) - J*K';
S11 =  G(:,:,1) - K*K';

% Eigenvalues and eigenvectors

W = S01'/chol(S00);
[V,D] = eig(W*W',S11,'chol');
[lam,sidx] = sort(diag(D),'descend');
assert(all(lam>=0 & lam<1),'Bad eigenvalues!');
V = V(:,sidx); % eigenvectors
if normevs
	VS11 = chol(S11)*V;
	V = V./sqrt(diag(VS11'*VS11))';
end
loglam = log(1-lam);

% Statistical results

stats.evs =  lam;
stats.A   =  S01*V;
stats.B   =  V;
stats.me  = -loglam;                         % scale by effective sample size for inference
stats.tr  = -flipud(cumsum(flipud(loglam))); % scale by effective sample size for inference
