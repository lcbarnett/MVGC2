% SINUSOIDAL LEAST-SQUARES-FIT FREQUENCY
%
% INVOCATION
%
% ffit = sinufit_ls(x,fs,f,fftol,verb)
%
% mse  = sinufit_ls(x,fs,f,'MSE')
%
% INPUTS
%
% x         matrix of time-series values (variables x observations)
% fs        sampling frequency (default: angular frequency in range 0 ... 2*pi)
% f         for the first form, an ascending 2-vector specifying the fit frequency range; for
%           the second 'MSE' form, a vector of frequencies at which the MSE will be calculated.
% fftol     frequency fit tolerance (default: 1e-8)
% verb      verbosity
%
% OUTPUTS
%
% ffit      fitted frequency
% mse       the MSE at the supplied frequencies
%
% NOTE      If calculating the optimal frequency in a range, the range should be
%           small enough that the MSE is approximately "U-shaped" within that
%           range. If not, the 'fminbnd' function used to locate the optimal
%           frequency may return a value corresponding to a local sub-optimum.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function retval = sinufit_ls(x,fs,f,fftol,verb)

if isempty(fs), fs = 2*pi; end % If no sampling frequency supplied assume angular frequency.

cmse = nargin > 3 && ischar(fftol) && strcmpi(fftol,'MSE');

assert(ismatrix(x),'Time series must be a matrix');
[n,m] = size(x); % n = number of variables, m = number of observations

if cmse
	assert(nargin == 4,'Too many input arguments for ''MSE'' option');
	assert(isvector(f),'For ''MSE'' option, a vector of frequencies must be supplied');
	f = f(:); % ensure column vector
else
	assert(nargin > 2,'Too few input arguments');
	if nargin < 4 || isempty(fftol), fftol = 1e-8;  end % Default tolerance for 'fminbnd' frequency fit
	if nargin < 5 || isempty(verb),  verb  = false; end % Verbosity (if set, print 'fminbnd' diagnostics)
	assert(isvector(f) && length(f) == 2 && f(2) >= f(1),'Fit frequencies must be an ascending 2-vector');
end
nf = size(f,1);

mu  = mean(x,2);     % save mean
sig = std(x,[],2);   % save std. dev.
x   = (x-mu)./sig;   % demean and normalise by variance
t   = 0:(m-1);       % time sequence
w   = (2*pi)*(f/fs); % convert to angular frequencies

if cmse
	nf = length(f);
	retval = zeros(nf,n);
	for k = 1:nf
		retval(k,:) = MSES(w(k));
	end
else
	if verb
		os = optimset('Display','iter','TolX',fftol);
	else
		os = optimset('TolX',fftol);
	end
	[wfit,~,exitflag] = fminbnd(@MSE,w(1),w(2),os); % find wfit which minimises the MSE

	if exitflag == 1,
		retval = (fs*wfit)/(2*pi); % fit frequency (convert back to ordinary frequency)
	else
		retval = NaN;
	end
end

% Nested functions to calculate MSE

    function mse = MSE(ww)

		% Passed as parameter to 'fminbnd' (see above). The MSE is calculated up to
		% additive/multiplicative factors which don't depend on w.

		s = [sin(ww*t);cos(ww*t)]; % sinusoids
		e = x-((x*s')/(s*s'))*s;   % error (OLS)
		mse = mean(mean(e.^2,2));  % mean of mean square errors across variables

    end % function MSE

    function mses = MSES(ww)

		% MSEs per variable

		s = [sin(ww*t);cos(ww*t)]; % sinusoids
		e = x-((x*s')/(s*s'))*s;   % error (OLS)
		mses = mean(e.^2,2);       % mean square errors for all variables

    end % function MSES

end % function sinufit
