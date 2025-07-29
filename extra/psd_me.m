function [S,f] = psd_me(X,fs,leakage)

% Power spectral density for multi-epoched data

if nargin < 2, fs      = []; end % pspectrum default 2*pi
if nargin < 3, leakage = []; end % pspectrum default 0.5

X = permute(X,[2,3,1]); % obs x epochs x chans

for i = 1:size(X,3)
	[Si,f] = pspectrum(X(:,:,i),fs,'power','Leakage',leakage);
	S(:,i) = mean(Si,2); % average over epochs
end
