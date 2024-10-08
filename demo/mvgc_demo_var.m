%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test data generation

if ~exist('ntrials',   'var'), ntrials   = 5;       end % number of trials (1 for single-trial)
if ~exist('nobs',      'var'), nobs      = 500;     end % number of observations per trial
if ~exist('fs',        'var'), fs        = 200;     end % sample rate (Hz)

% VAR model parameters

if ~exist('tnet',      'var'), tnet      = tnet5;   end % connectivity graph
if ~exist('moact',     'var'), moact     = 6;       end % model order
if ~exist('rho',       'var'), rho       = 0.95;    end % spectral radius
if ~exist('wvar',      'var'), wvar      = 0.5;     end % var coefficients decay weighting factor
if ~exist('rmi',       'var'), rmi       = 0.5;     end % residuals log-generalised correlation (multi-information):
                                                        % rmi = -log|R|. rmi = 0 yields zero correlation, rmi = [] is
                                                        % uniform random on space of correlation matrices
% VAR model order estimation

if ~exist('moregmode', 'var'), moregmode = 'LWR';   end % VAR model estimation regression mode ('OLS' or 'LWR')
if ~exist('mosel',     'var'), mosel     = 'LRT';   end % model order selection ('ACT', 'AIC', 'BIC', 'HQC', 'LRT', or supplied numerical value)
if ~exist('momax',     'var'), momax     = 2*moact; end % maximum model order for model order selection

% VAR model parameter estimation

if ~exist('regmode',   'var'), regmode   = 'LWR';   end % VAR model estimation regression mode ('OLS' or 'LWR')

% MVGC (time domain) statistical inference

if ~exist('alpha',     'var'), alpha     = 0.05;    end % significance level for Granger casuality significance test
if ~exist('stest',     'var'), stest     = 'F';     end % statistical inference test: 'F' or 'LR' (likelihood-ratio chi^2)
if ~exist('mhtc',      'var'), mhtc      = 'FDRD';  end % multiple hypothesis test correction (see routine 'significance')

% MVGC (frequency domain)

if ~exist('fres',      'var'), fres      = [];      end % spectral MVGC frequency resolution (empty for automatic calculation)

%%%%%%%%%%%%%%%%%%%%%%%%% Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('seed',      'var'), seed      = 0;       end % random seed (0 for unseeded)
if ~exist('plotm',     'var'), plotm     = 0;       end % plot mode (figure number offset, or Gnuplot terminal string)

%%%%%%%%%%%%%%%%%%%%%%%%% Test data generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate random VAR test data
%
% NOTE: This is where you would read in your own time series data; it should
% be assigned to the 3D array X (variables x observations x trials) if multi-
% trial, or 2D matrix (variables x observations) if single-trial.

% Seed random number generator.

rng_seed(seed);

% Generate random VAR coefficients for test network.

AA = var_rand(tnet,moact,rho,wvar);
nvars = size(AA,1); % number of variables

% Generate random residuals covariance (in fact correlation) matrix.

VV = corr_rand(nvars,rmi);

% Report information on the generated VAR model and check for errors.

infoo = var_info(AA,VV);
assert(~infoo.error,'VAR error(s) found - bailing out');

% Generate multi-trial VAR time series data with normally distributed residuals
% for generated VAR coefficients and residuals covariance matrix.

ptic('*** var_to_tsdata... ');
X = var_to_tsdata(AA,VV,nobs,ntrials,'stationary');
ptoc;

%%%%%%%%%%%%%%%%%%%%%%%%% VAR modelling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove temporal mean and normalise by temporal variance.
% Not strictly necessary, but may help numerical stability
% if data has very large or very small values.

X = demean(X,true);

% Model order estimation

% Calculate and plot VAR model order estimation criteria up to specified maximum model order.

ptic('\n*** tsdata_to_varmo... ');
if isnumeric(plotm), plotm = plotm+1; end
[moaic,mobic,mohqc,molrt] = tsdata_to_varmo(X,momax,moregmode,[],[],plotm);
ptoc;

% Select and report VAR model order.

morder = moselect(sprintf('VAR model order selection (max = %d)',momax), mosel,'ACT',moact,'AIC',moaic,'BIC',mobic,'HQC',mohqc,'LRT',molrt);
assert(morder > 0,'selected zero model order! GCs will all be zero!');
if morder >= momax, fprintf(2,'*** WARNING: selected maximum model order (may have been set too low)\n'); end

% VAR model estimation

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,V] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed - bailing out');

% Report information on the estimated VAR, and check for errors.
%
% IMPORTANT: We check the VAR model for stability and symmetric positive-
% definite residuals covariance matrix. THIS CHECK IS ESSENTIAL - subsequent
% routines may fail if there are errors here. If there are problems with the
% data (e.g. non-stationarity, colinearity, etc.) there's also a good chance
% they'll show up at this point - and the diagnostics may supply useful
% information as to what went wrong.

info = var_info(A,V);
assert(~info.error,'VAR error(s) found - bailing out');

%%%%%%%%%%%%%%%%%%% Granger causality analysis - time domain %%%%%%%%%%%%%%%%%%%

% Estimate time-domain pairwise-conditional GC (information transfer rate in nats/second)

ptic('*** var_to_pwcgc... ');
F = var_to_pwcgc(A,V);
ptoc;
assert(~isbad(F,false),'GC estimation failed');

% Calculate test statistics (F or likelihood-ratio chi^2 test) for time-domain pairwise-conditional GC hypothesis test

ptic('*** var_to_pwcgc_tstat... ');
tstat = var_to_pwcgc_tstat(X,V,morder,regmode,stest);
ptoc;

% Calculate p-values for test statistics

pval = mvgc_pval(tstat,stest,1,1,nvars-2,morder,nobs,ntrials); % for pairwise-conditional, nx = 1, ny = 1, nz = nvars-2

% Significance test p-values (F or likelihood-ratio chi^2 test), correcting for multiple hypotheses.

sig = significance(pval,alpha,mhtc);

% For comparison, we also calculate the actual pairwise-conditional causalities

ptic('*** var_to_pwcgc... ');
FF = var_to_pwcgc(AA,VV);
ptoc;
assert(~isbad(FF,false),'GC calculation failed');

% Plot time-domain causal graph and significance.

maxF = 1.1*max(nanmax(F(:),nanmax(FF(:))));
pdata = {FF,F,sig};
ptitle = {'PWCGC (actual)','PWCGC (estimated)',[stest '-test']};
maxp = [maxF,maxF,1];
if isnumeric(plotm), plotm = plotm+1; end
plot_gc(pdata,ptitle,[],maxp,plotm,[0.6,2.5]);

%%%%%%%%%%%%%%%% Granger causality analysis - frequency domain %%%%%%%%%%%%%%%%%

% If not specified, we set the frequency resolution to something sensible.

if isempty(fres)
	ptic('\n*** var2fres... ');
	[fres,frierr,frpow2] = var2fres(A,V);
	ptoc;
	fprintf('\nUsing frequency resolution %d = 2^%d (integration error = %.2e)\n',fres,frpow2,frierr);
else
	frierr = var_check_fres(A,V,fres);
	fprintf('\nUsing frequency resolution %d (integration error = %.2e)\n',fres,frierr);
end

ptic(sprintf('\n*** var_to_spwcgc (at frequency resolution = %d)... ',fres));
f = var_to_spwcgc(A,V,fres);
ptoc;
assert(~isbad(f,false),'spectral GC estimation failed');

% For comparison, we also calculate the actual pairwise-conditional spectral causalities

ptic(sprintf('*** var_to_spwcgc (at frequency resolution = %d)... ',fres));
ff = var_to_spwcgc(AA,VV,fres);
ptoc;
assert(~isbad(ff,false),'spectral GC calculation failed');

% Get frequency vector corresponding to specified sampling rate.

freqs = sfreqs(fres,fs);

% Plot spectral causal graphs.

if isnumeric(plotm), plotm = plotm+1; end
plot_sgc({ff,f},freqs,'Spectral Granger causalities (blue = actual, red = estimated)',[],plotm);

% Granger causality calculation: frequency domain -> time-domain

% Check that spectral causalities average (integrate) to time-domain
% causalities. Note that this may occasionally fail if a certain condition
% on the VAR parameters is not satisfied (Geweke 1982). This does not
% invalidate the time-domain GC!

Fint = bandlimit(f,3); % integrate spectral MVGCs (frequency is dimension 3 of CPSD array)

fprintf('\n*** GC spectral integral check... ');
rr = abs(F-Fint)./(1+abs(F)+abs(Fint)); % relative residuals
mrr = max(rr(:));                       % maximum relative residual
if mrr < 1e-6
    fprintf('PASS: max relative residual = %.2e\n',mrr);
else
    fprintf(2,'FAIL: max relative residual = %.2e (too big!)\n',mrr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
