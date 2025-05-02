%%%%%%%%%%%%%%%%%% Default parameters - override on command line %%%%%%%%%%%%%%%

% Model order estimation

if ~exist('empdata', 'var'), empdata  = false; end % set this to true if you have empirical time-series data (see below)
if ~exist('varmosel','var'), varmosel = 'HQC'; end % VAR model order selection ('AIC', 'BIC', 'HQC', 'LRT', or supplied numerical value)
if ~exist('varmomax','var'), varmomax = 48;    end % maximum model order for VAR model order selection
if ~exist('ernorm',  'var'), ernorm   = false; end % calculate normalised entropy rates (this makes them scale-invariant - and always negative!)
if ~exist('fres',    'var'), fres     = [];    end % frequency resolution for spectral entropy rates; leave empty for automatic calculation

%%%%%%%%%%%%%%%%%% Read your data in our generate test data %%%%%%%%%%%%%%%%%%%%

if empdata % You have empirical time-series data

	% Your data should be read into a 3d array called 'X' before running this script
	% (or a 2d matrix if it is not epoched). The dimensions of X must be:
	%
	%    (number of channels) x (number of observations) x (number of epochs)
	%
	% or just
	%
	%    (number of channels) x (number of observations)
	%
	% if not epoched. The data should be in double-precision floating-point format.
	% You will also need to set the variable 'fs' to the sampling rate for your data.

	[nchans,nobs,nepochs] = size(X); % this is still okay if X is 2d (not epoched)

	% fs = ???; % set appropriately

else  % Generate test data from a random state-space model

	% State-space test data generation parameters

	if ~exist('nchans',  'var'), nchans   = 5;     end % number of channels
	if ~exist('nobs',    'var'), nobs     = 2000;  end % number of observations per trial
	if ~exist('nepochs', 'var'), nepochs  = 10;    end % number of trials
	if ~exist('ssmo',    'var'), ssmo     = 10;    end % SS model order
	if ~exist('rhoa',    'var'), rhoa     = 0.9;   end % AR spectral radius
	if ~exist('rmi',     'var'), rmi      = 0.5;   end % residuals log-generalised correlation (multi-information)
	if ~exist('fs',      'var'), fs       = 200;   end % sample rate (Hz)
	if ~exist('seed',    'var'), seed     = 0;     end % random seed for reproducibility (0 for unseeded)

	% Seed random number generator.

	rng_seed(seed);

	% Generate random SS parameters in innovations form

	[AA,CC,KK] = iss_rand(nchans,ssmo,rhoa);

	% Generate random residuals covariance (correlation) matrix.

	VV = corr_rand(nchans,rmi);

	% Report information on the generated SS model and check for errors.

	infoo = ss_info(AA,CC,KK,VV);
	assert(~infoo.error,'SS error(s) found - bailing out');

	% Generate multi-trial SS time series data with normally distributed residuals
	% for ISS parameters (AA,CC,KK,VV)

	ptic('*** ss_to_tsdata... ');
	X = ss_to_tsdata(AA,CC,KK,VV,nobs,nepochs);
	ptoc;

end

%%%%%%%%%%%%%%% Estimate SS model from data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate and plot VAR model order estimation criteria up to specified maximum model order.
% NOTE: VAR model order is required for SS-SS estimation.

ptic('\n*** tsdata_to_varmo... ');
fignum = 1;
[varmoaic,varmobic,varmohqc,varmolrt] = tsdata_to_varmo(X,varmomax,'LWR',[],[],fignum);
ptoc;

% Select and report VAR model order. NOTE: may sometimes fail - this is normal :-/

varmo = moselect(sprintf('VAR model order selection (max = %d)',varmomax),varmosel,'AIC',varmoaic,'BIC',varmobic,'HQC',varmohqc,'LRT',varmolrt);
assert(varmo > 0,'selected zero model order!');
if varmo >= varmomax, fprintf(2,'*** WARNING: selected VAR maximum model order (may have been set too low)\n'); end

% Estimate SS model using CCA SS-SS algorithm with SVC model order criterion

ptic('\n*** tsdata_to_ss... ');
fignum = fignum+1;
[A,C,K,V] = tsdata_to_ss(X,2*varmo,'SVC',fignum); % Bauer recommends 2 x VAR model order for past/future horizon
ptoc;

% If the ernorm flag is set to true, then the ISS model is normalised by
% process covariance. In this case, error rates are equivalent to minus the
% mutual information between the process and its own past, and as such are
% scale-invariant (and non-positive).

if ernorm
	fprintf('\n*** Normalising ISS model by variance\n');
	[A,C,K,V] = ss_normalise(A,C,K,V);
end

% Report information on the estimated SS, and check for errors.

info = ss_info(A,C,K,V);
assert(~info.error,'SS model error(s) found - bailing out');

%%%%%%%%%%%%%%% Frequency resolution for spectral entropy rates %%%%%%%%%%%%%%%%

if isempty(fres)
	ptic('*** ss2fres... ');
	[fres,frierr,frpow2] = ss2fres(A,C,K,V);
	ptoc;
	fprintf('\n*** Using automatically-calculated frequency resolution %d = 2^%d (integration error = %.2e)\n',fres,frpow2,frierr);
else
	assert(isscalar(fres) && isnumeric(fres) && fres > 0,'Frequency resolution must be a positive number');
	fprintf('*** Using user-supplied frequency resolution %d\n',fres);
end

%%%%%%%%%%%%%%% Calculate entropy rates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time domain, global

ptic('\n*** ss_to_erate: time domain, global... ');
glerate = ss_to_erate(A,C,K,V,'allchans');
ptoc;

% Time domain, per-channel

ptic('\n*** ss_to_erate: time domain, per-channel... ');
pcerate = ss_to_erate(A,C,K,V,'perchan');
ptoc;

% Broadband, global

ptic('\n*** ss_to_serate: frequency domain (broadband), global... ');
glbrerate = ss_to_serate(A,C,K,V,'allchans','broadband',fs,fres);
ptoc;

% Broadband, per-channel

ptic('\n*** ss_to_serate: frequency domain (broadband), per-channel... ');
pcbrerate = ss_to_serate(A,C,K,V,'perchan','broadband',fs,fres);
ptoc;

% Standard frequency bands, global

ptic('\n*** ss_to_serate: frequency domain (standard bands), global... ');
glsterate = ss_to_serate(A,C,K,V,'allchans','stdx',fs,fres);
ptoc;

% Standard frequency bands, per-channel

ptic('\n*** ss_to_serate: frequency domain (standard bands), per-channel... ');
pcsterate = ss_to_serate(A,C,K,V,'perchan','stdx',fs,fres);
ptoc;

% Display time-domain and frequency band-limited entropy rates in a table

fprintf('\n------------------------');           for i = 1:nchans, fprintf('-----------');               end; fprintf('\n');
fprintf('Entropy rates :   GLOBAL');             for i = 1:nchans, fprintf('     chan %d',i);            end; fprintf('\n');
fprintf('------------------------');             for i = 1:nchans, fprintf('-----------');               end; fprintf('\n');
fprintf('TIME DOMAIN   :  % 7.4f',glerate);      for i = 1:nchans, fprintf('    % 7.4f',pcerate(i));     end; fprintf('\n');
fprintf('delta         :  % 7.4f',glbrerate(1)); for i = 1:nchans, fprintf('    % 7.4f',pcsterate(1,i)); end; fprintf('\n');
fprintf('theta         :  % 7.4f',glbrerate(2)); for i = 1:nchans, fprintf('    % 7.4f',pcsterate(2,i)); end; fprintf('\n');
fprintf('alpha         :  % 7.4f',glbrerate(3)); for i = 1:nchans, fprintf('    % 7.4f',pcsterate(3,i)); end; fprintf('\n');
fprintf('beta          :  % 7.4f',glbrerate(4)); for i = 1:nchans, fprintf('    % 7.4f',pcsterate(4,i)); end; fprintf('\n');
fprintf('low-gamma     :  % 7.4f',glbrerate(5)); for i = 1:nchans, fprintf('    % 7.4f',pcsterate(5,i)); end; fprintf('\n');
fprintf('high-gamma    :  % 7.4f',glbrerate(6)); for i = 1:nchans, fprintf('    % 7.4f',pcsterate(6,i)); end; fprintf('\n');
fprintf('------------------------');             for i = 1:nchans, fprintf('-----------');               end; fprintf('\n');

% Sanity check: broadband/band-limited entropy rates (global and per-channel) should
% integrate/sum respectively to the correpsonding time-domain entropy rates.

glbrinterror = abs(bandlimit(glbrerate)-glerate);
pcbrinterror = maxabs(bandlimit(pcbrerate)-pcerate);
glstinterror = max(abs(sum(glsterate)-glerate));
pcstinterror = max(abs(sum(pcsterate)-pcerate'));
maxabserror  = max([glbrinterror,pcbrinterror,glstinterror,pcstinterror]);
fprintf('\nSANITY CHECK: max. absolute integration error = %.1e\n\n',maxabserror);

% Plot broadband spectral entropy rates

freqs = (fs/2)*(0:fres)'/fres; % vector of frequencies up to Nyqvist

fignum = fignum+1;
figure(fignum); clf;
plot(freqs,glbrerate,'k','linewidth',1.5);
hold on
plot(freqs,pcbrerate);
title('Broadband spectral entropy rates');
xlabel('Frequency (Hz)');
ylabel('Entropy rate');
legs = cell(nchans+1,1);
legs{1} = 'GLOBAL';
for i = 1:nchans,legs{i+1} = sprintf('channel %2d',i); end
legend(legs,'location','northeastoutside');
hold off
grid on
