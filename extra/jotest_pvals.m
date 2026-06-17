function pval = jotest_pvals(n,cvtab,sigvals,tstat)

% P-values values for Johansen tests with n variables for test
% statistics tstat, for cointegration rank 0 .. n-1.
%
% Call
%
%	[cvtab,sigvals] = jotest_cvtab(model)
%
% to get appropriate critical value tables and significance levels.
%
% NOTE: critical values only available for last 12 cointegration ranks;
% lower values returned as NaNs.

rmin = max(n-12,0);

pval.tr = nan(n,1);
for r = rmin:n-1
	if     tstat.tr(r+1) >= cvtab(n-r,1,  1), pval(r+1) = sigvals(1);
	elseif tstat.tr(r+1) <= cvtab(n-r,end,1), pval(r+1) = sigvals(end);
	else   pval.tr(r+1) = interp1(cvtab(n-r,:,1),sigvals,tstat.tr(r+1),'linear');
	end
end

pval.me = nan(n,1);
for r = rmin:n-1
	if     tstat.me(r+1) >= cvtab(n-r,1,  2), pval(r+1) = sigvals(1);
	elseif tstat.me(r+1) <= cvtab(n-r,end,2), pval(r+1) = sigvals(end);
	else   pval.me(r+1) = interp1(cvtab(n-r,:,2),sigvals,tstat.me(r+1),'linear');
	end
end
