function stats = jotest_tsdata(y,p,normevs,lamtol)

% Calculate VECM Johansen Test statistics from multitrial time-series data
%
% This is basically the same as the Matlab 'jcitest' (econ. Toolbox), but
% accommodates multi-trial data, and works for n > 12.
%
% NOTE: p is the VAR, not VECM autoregressive order!

if nargin < 3 || isempty(normevs), normevs = true;  end % normalise eigenvectors?
if nargin < 4 || isempty(lamtol),  lamtol  = 1e-10; end % eigenvalue tolerance

[n,T,N] = size(y);

DY  = y(:,2:T,:)-y(:,1:T-1,:); % differenced time series
LY  = y(:,1:T-1,:);            % lagged time series
T   = T-1;                     % differencing/lagging loses 1st observation!
q   = p-1;                     % VECM autoregressive order
obs = p:T;                     % useable observations (i.e., lose the first q)
Te  = N*(T-q);                 % effective number of observations (sample size)

% Stack q lags of differences (these are the regressors)

LDY = zeros(n,q,Te); % lagged differences
for k = 1:q
	LDY(:,k,:) = reshape(DY(:,obs-k,:),n,Te); % concatenate trials for k-lagged differences
end
LDY = reshape(LDY,q*n,Te); % stack lagged differenced time series

% Align DY, LY with LDY

DY = reshape(DY(:,obs,:),n,Te); % concatenate trials for differenced time series
LY = reshape(LY(:,obs,:),n,Te); % concatenate trials for lagged time series

% Residuals of DY and LY regressions on LDY

R0 = DY - (DY/LDY)*LDY;
R1 = LY - (LY/LDY)*LDY;

% (Cross-)covariance matrices

S00 = (R0*R0')/Te;
S01 = (R0*R1')/Te;
S11 = (R1*R1')/Te;

% Eigenvalues and eigenvectors

W = S01'/chol(S00);
[V,D] = eig(W*W',S11,'chol');
[lam,sidx] = sort(diag(D),'descend');
assert(all(lam>=0 & lam < 1+lamtol),'Bad eigenvalues!');
lam(lam > 1) = 1;
V = V(:,sidx); % eigenvectors
if normevs
	VS11 = chol(S11)*V;
	V = V./sqrt(diag(VS11'*VS11))';
end
loglam = log(1-lam);

% Statistical results

stats.ess =  Te; % effective sample size
stats.evs =  lam;
stats.A   =  S01*V;
stats.B   =  V;
stats.me  = -loglam;                         % scale by effective sample size for inference
stats.tr  = -flipud(cumsum(flipud(loglam))); % scale by effective sample size for inference
