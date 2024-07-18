function [S,f,window,overlap,nfft] = cpsd_me(X,fs,window,overlap,nfft,autospec,verb)

% Welch method CPSD or PSD for multi-epoched data

[nchans,nobs,nepochs] = size(X);

if nargin < 3 || isempty(window),   window   = round(nobs/2);    end
if nargin < 4 || isempty(overlap),  overlap  = round(window/2);  end
if nargin < 5 || isempty(nfft),     nfft     = 2^nextpow2(nobs); else, assert(mod(nfft,2)==0,'nfft must be even!'); end
if nargin < 6 || isempty(autospec), autospec = false;            end
if nargin < 7 || isempty(verb),     verb     = true;             end

f = (((0:nfft/2)')/nfft)*fs; % frequencies on [0,fs/2]

hfft = nfft/2+1;
if autospec
	S = zeros(hfft,nchans);
else
	S = zeros(hfft,nchans,nchans);
end

X = permute(X,[2 1 3]);

if autospec % NOTE: Autospectra (PSDs) matrix S is frequencies x channels
	for e = 1:nepochs
		S = S + pwelch(X(:,:,e),window,overlap,nfft,fs);
	end
	S = S/nepochs;
	S(1  ,:) = 2*S(1,  :);
	S(end,:) = 2*S(end,:);

else        % NOTE: Cross-spectra (CPSDs) array S is channels x channels x frequencies
	for e = 1:nepochs
		if verb, fprintf('\tepoch %3d of %3d\n',e,nepochs); end
		for i = 1:nchans
			S(:,:,i) = S(:,:,i) + cpsd(X(:,i,e),X(:,:,e),window,overlap,nfft,fs);
		end
	end
	S = permute(S/nepochs,[3 2 1]);
	S(:,:,1  ) = 2*S(:,:,1  );
	S(:,:,end) = 2*S(:,:,end);
end
