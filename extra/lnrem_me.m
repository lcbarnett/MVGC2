%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Remove line noise from multi-channel, multi-epoch data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Line noise frequencies are firstly fitted per channel. An
% appropriate PSD at the line noise frequency is then interpolated
% by a least-squares log-log fit around the line noise frequency,
% leaving a small gap around the line noise frequency. For each channel,
% sinusoids are then fitted for each epoch using least-squares; the
% weighted sinusoids are then subtracted. Suitable per-channel weightings
% are calculated via a binary chop, so that the channel PSD matches the
% interpolated value at the line noise frequency.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = lnrem_me( ...
	X,       ... % Time series data: channels x observations x epochs
	fs,      ... % sampling rate (Hz)
	lnfreq,  ... % line noise frequency (Hz)
	nharms,  ... % number of line noise harmonics
	welch,   ... % Use Welch method (else periodogram)
	window,  ... % Welch spectral estimation window size (observations)
	overlap, ... % Welch spectral estimation overlap (observations)
	nfft,    ... % Number of FFT data points
	ffitr,   ... % frequency fit range (Hz)
	fftol,   ... % frequency fit tolerance
	finti,   ... % frequency interpolation interval (Hz)
	fintig,  ... % frequency interpolation interval gap (Hz)
	bctol,   ... % binary chop tolerance
	bcmaxi   ... % binary chop maximum iterations
)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nchans,nobs,nepochs] = size(X);

% Parameter defaults (X, fs, lnfreq and nharms are mandatory)

if nargin <  4 || isempty(nharms),  nharms  = floor(fs/2/lnfreq); end
if nargin <  5 || isempty(welch),   welch   = true;               end
if nargin <  6 || isempty(window),  window  = round(nobs/2);      end
if nargin <  7 || isempty(overlap), overlap = round(window/2);    end
if nargin <  8 || isempty(nfft),    nfft    = 2^nextpow2(nobs);   end
if nargin <  9 || isempty(ffitr),   ffitr   = 0.5;                end
if nargin < 10 || isempty(fftol),   fftol   = 1e-8;               end
if nargin < 11 || isempty(finti),   finti   = 20;                 end
if nargin < 12 || isempty(fintig),  fintig  = 1.5;                end
if nargin < 13 || isempty(bctol),   bctol   = 1e-10;              end
if nargin < 14 || isempty(bcmaxi),  bcmaxi  = 100;                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nlnrem_me: removing line noise\n');

% Frequency and time vectors

hfft  = nfft/2+1;
f     = ((1:hfft)'/hfft)*(fs/2); % zero to Nyqvist  frequencies
logf  = log(f);
t     = (0:nobs-1)'; % time stamps

for hnum = 1:nharms % loop through harmonics
	lnf = hnum*lnfreq;

	% Fit line noise frequency per epoch (i.e., across channels) using OLS;
	% the logic is that at any given epoch, the line noise frequency should
	% be the same for all channels (although phase and magnitude may vary
	% across channels)

	fprintf('\nHarmonic %d of %d = %5.2fHz : fitting per-epoch frequencies\n',hnum,nharms,lnf);
	ffits = zeros(nepochs,1);
	for e = 1:nepochs % loop through epochs
		ffits(e) = sinufit_ls(X(:,:,e),fs,[lnf-ffitr,lnf+ffitr],fftol); % in MVGC2/extra
		fprintf('\tepoch %3d of %3d : lnfreq = %6.2fHz\n',e,nepochs,ffits(e));
	end % epochs

	% Regression local frequency range parameters (leave gap around line noise frequency)

	iffit = findfreq(f,lnf); % index of nearest frequency to fitted line noise frequency in frequency vector f
	ireg  = (f >= lnf-finti & f <= lnf-fintig) | (f <= lnf+finti & f >= lnf+fintig); % indices of local range frequencies
	ivar  = [logf(ireg) ones(nnz(ireg),1)]; % log-log fit regressor: log-frequency in local range

	fprintf('\nHarmonic %d of %d = %5.2fHz  : removing line-noise\n',hnum,nharms,lnf);
	for i = 1:nchans % loop through channels

		fprintf('\tchannel %2d of %2d : ',i,nchans);

		% De-mean and normalise (save means and std. devs.)

		x   = squeeze(X(i,:,:)); % observations x epochs
		mu  = mean(x);
		sig = std(x);
		x   = (x-mu)./sig;

		% Calculate desired channel PSD (interpolate PSD locally by log-log fit)

		if welch
			psd = mean(pwelch(x,window,overlap,nfft,fs),2);
		else
			psd = mean(periodogram(x,[],nfft,'psd',fs),2);
		end
		lpsd    = log(psd);                  % log PSD
		dvar    = lpsd(ireg);                % log-log fit regresee: log PSD in local range
		regc    = (ivar'*ivar)\(ivar'*dvar); % regression coefficients (OLS)
		lpsdreg = regc(1)*logf+regc(2);      % regression line (local fit)
		psdreg  = exp(lpsdreg);              % regression PSD
		psdfit  = psdreg(iffit);             % channel PSD at line-noise freqency we want (target)

		% Calculate sinusoidal components per epoch, using per-epoch line noise frequencies

		s = zeros(nobs,nepochs);
		for e = 1:nepochs % loop through epochs
			ffit      = ffits(e);
			wfit      = ((2*pi)/fs)*ffit;           % circular frequency
			ss        = [sin(wfit*t) cos(wfit*t)];  % sinusoid
			s(:,e)    = ss*((ss'*ss)\(ss'*x(:,e))); % regress signal on sinusoid (OLS)
		end

		% Use binary chop to calculate sinusoid weighting to achieve interpolated
		% channel PSD, and subtract weighted sinusoids from signal across epochs

		wmax = 1;
		wmin = 0;
		w = (wmin+wmax)/2; % weighting for this channel: 0 <= w <= 1
		dold = Inf;
		k = 0;
		success = true;
		while true
			k = k+1;
			if k > bcmaxi, success = false; break; end % binary chop timed out (shouldn't happen!)
			y = x-w*s; % channel data with weighted sinusoids subtracted
			if welch
				ypsd = mean(pwelch(y,window,overlap,nfft,fs),2);
			else
				ypsd = mean(periodogram(y,[],nfft,'psd',fs),2);
			end
			ypsdfit = ypsd(iffit);
			d = abs(ypsdfit-psdfit);
			dd = abs(d-dold);
			dold = d;
			if d < bctol | dd < bctol, break; end % success: we have correct weighting
			if     ypsdfit < psdfit % lnf power too low: decrease weighting
				wmax = w;
			elseif ypsdfit > psdfit % lnf power too high: increase weighting
				wmin = w;
			end
			w = (wmin+wmax)/2;
		end

		fprintf('iters = %2d (%.2e, %.2e), weight = %6.4f\n',k,d,dd,w);

		X(i,:,:) = mu+sig.*y; % restore means and std. devs.

	end % channel

end % harmonics

function i = findfreq(f,freq) % Return index of value nearest to freq in frequency vector f

	i1 = find(f <= freq,1,'last' );
	i2 = find(f >= freq,1,'first');
	if 2*freq <= f(i1)+f(i2)
		i = i1;
	else
		i = i2;
	end

end
