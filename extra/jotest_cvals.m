function cval = jotest_cvals(n,cvtab,sigvals,alpha)

% Critical values for Johansen tests with n variables at significance
% level alpha, for cointegration rank 0 .. n-1.
%
% Call
%
%	[cvtab,sigvals] = jotest_cvtab(model)
%
% to get appropriate critical value tables and significance levels.
%
% NOTE 1: critical values only available for last 12 cointegration ranks;
% lower values returned as NaNs.
%
% NOTE 2: to test the appropriate null hypotheses, the test statistics
% (calculated by jotest_tsdata) must be scaled by effective sample size.

if nargin < 4 || isempty(alpha), alpha = 0.05; end % significance level

rmin = max(n-12,0);

cval.tr = nan(n,1);
for r = rmin:n-1
	cval.tr(r+1) = interp1(sigvals,cvtab(n-r,:,1),alpha,'linear');
end

cval.me = nan(n,1);
for r = rmin:n-1
	cval.me(r+1) = interp1(sigvals,cvtab(n-r,:,2),alpha,'linear');
end
