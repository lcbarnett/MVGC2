function cval = jotest_cvals(n,model,alpha)

% Critical values for Johansen tests for given model with n variables
% at significance level alpha, for cointegration rank 0 .. n-1.
%
% NOTE 1: critical values only available for last 12 cointegration ranks;
% lower values returned as NaNs.
%
% NOTE 2: to test the appropriate null hypotheses, the test statistics
% (calculated by jotest_tsdata) must be scaled by effective sample size.

if nargin < 3 || isempty(alpha), alpha = 0.05; end % significance level

CV = load('Data_JCITest');
switch model
	case 'H2',  cvtab = squeeze(CV.JCV(:,:,1,:));
	case 'H1*', cvtab = squeeze(CV.JCV(:,:,2,:));
	case 'H1',  cvtab = squeeze(CV.JCV(:,:,3,:));
	case 'H*',  cvtab = squeeze(CV.JCV(:,:,4,:));
	case 'H',   cvtab = squeeze(CV.JCV(:,:,5,:));
end

sigvals = [0.001 0.005:0.005:0.10 0.125:0.025:0.20 0.80:0.025:0.875 0.90:0.005:0.995 0.999];

rmin = max(n-12,0);

cval.tr = nan(n,1);
for r = rmin:n-1
	cval.tr(r+1) = interp1(sigvals,cvtab(n-r,:,1),alpha,'linear');
end

cval.me = nan(n,1);
for r = rmin:n-1
	cval.me(r+1) = interp1(sigvals,cvtab(n-r,:,2),alpha,'linear');
end
