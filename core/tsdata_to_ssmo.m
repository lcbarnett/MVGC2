function [mosvc,rmax] = tsdata_to_ssmo(X,pf,plotm)

% Estimate SS model porder using Bauer's Singular Value Criterion (SVC; D. Bauer, Automatic 37, 2001)
%
% X       - observation process time series
% pf      - past/future horizons for canonical correlations
% plotm   - empty: don't plot, integer: Matlab plot to figure n (if zero, use next); string: Gnuplot terminal (may be empty)
%
% mosvc   - Bauer's Singular Value Criterion (SVC) optimal model order
% rmax    - maximum possible model order given supplied past/future horizons
%
% The past/future horizons pf may be supplied as a 2-vector [p,f] or a scalar p
% = f = pf. Bauer (D. Bauer, Automatic 37, 2001) recommends setting p = f = 2*p,
% where p is the optimal VAR model order for the observation process X according
% to Aikaike's Information Criterion (AIC).

[n,m,N] = size(X);

assert(all(isint(pf(:))),'past/future horizon must be a 2-vector or a scalar positive integer');
if isscalar(pf)
    p = pf;    f = pf;
elseif isvector(pf) && length(pf) == 2
    p = pf(1); f = pf(2);
else
    error('past/future horizon must be a 2-vector or a scalar positive integer');
end
assert(p+f < m,'past/future horizon too large (or not enough data)');
rmax = n*min(p,f);

X = demean(X); % no constant term (don't normalise!)

mp  = m-p;
mp1 = mp+1;
mf  = m-f;
mh  = mp1-f; % m-p-f+1

M  = N*mp;
M1 = N*mp1;
Mh = N*mh;

Xf = zeros(n,f,mh,N);
for k = 1:f
    Xf(:,k,:,:) = X(:,p+k:mf+k,:);
end
Xf = reshape(Xf,n*f,Mh);

XP = zeros(n,p,mp1,N);
for k = 0:p-1
    XP(:,k+1,:,:) = X(:,p-k:m-k,:);
end
Xp = reshape(XP(:,:,1:mh,:),n*p,Mh);
XP = reshape(XP,n*p,M1);

[Wf,cholp] = chol((Xf*Xf')/Mh,'lower');
assert(cholp == 0,'forward weight matrix not positive definite');

[Wp,cholp] = chol((Xp*Xp')/Mh,'lower');
assert(cholp == 0,'backward weight matrix not positive definite');

BETA = Xf/Xp; % 'OH' estimate: regress future on past
assert(all(isfinite(BETA(:))),'subspace regression failed');

[~,S] = svd(Wf\BETA*Wp); % SVD of CCA-weighted OH estimate

sval = diag(S);       % the singular values
df   = 2*n*(1:rmax)'; % number of free parameters (Hannan & Deistler, see also Bauer 2001) ... or rmax*rmax+2*n*rmax ???
svc  = -log(1-[sval(2:rmax);0]) + df*(log(Mh)/Mh); % Bauer's Singular Value Criterion

morder = (0:rmax)';
[~,idx] = min(svc); mosvc = morder(idx);

if ~isempty(plotm)
	plot_svc(sval,svc,mosvc,rmax,plotm);
end
