function erates = var_to_serate(A,V,rois,fbands,fs,fres)

% Calculate spectral entropy rates for SS model for specified
% frequency bands and ROIs.
%
% A is the VAR coefficients array, V the residuals covariance matrix
%
% rois is a cell vector, where each cell is a vector of channel
% numbers for the channels in an ROI. It may also take two special
% values: 'allchans' treats the entire system as a single roi;
% 'perchan' returns entropy rates for each individual channel.
%
% fbands is either a column vector of frequency band boundaries,
% a 2-column matrix of frequency bands, 'broadband', 'std' or 'stdx'.
% If set to "broadband", entropy rates across the entire spectrum)
% up to the Nyqvist frequency are returned. If 'std', then delta, theta,
% alpha, gamma frequency bands as defined in the neuroscience literature
% are used; 'stdx' splits gamma into high and low gamma.
%
% fs is sampling rate
%
% fres specifies the spectral resolution (number of frequencies,
% evenly spaced, at which spectral power is calculated). See
% ss2fres.m for automatic calculation of a frequency resolution for
% which spectral error rates will integrate/sum accurately to the
% corresponding time-domain error rate.
%
% Results are returned in the vector erates. In the case that fbands
% is specified as a vector of band boundaries, 'std', or 'stdx', as
% a reality check you can test whether sum(erates) == logdet(V) (if
% it's out by a couple of decimal points that's fine). Broadband
% entropy rates should also integrate to the corresponding logdet(V).

broadband = ischar(fbands) && strcmpi(fbands,'broadband');

if broadband
%	fprintf('\nBroadband: 0-%g Hz\n',fs/2);
else

	% Parse frequency bands specification

	if ischar(fbands)
		if     strcmpi(fbands,'std') % single gamma band > 30 Hz
			fbands = [4;8;15;30];
		elseif strcmpi(fbands,'stdx') % gamma band split into low (30 - 50 Hz) and high (> 50 Hz)
			fbands = [4;8;15;30;50];
		else
			error('Unknown frequency bands spec');
		end
	else
		assert(ismatrix(fbands) && (size(fbands,2) == 1 | size(fbands,2) == 2),'Frequency bands must be a string spec or a 1 or 2 column matrix');
	end
	if size(fbands,2) == 1 % vector of band boundaries
		assert(issorted(fbands),'Frequency boundaries must be sorted ascending');
		fb = fbands;
		nbb = size(fbands,1);
		fbands = zeros(nbb+1,2);
		fbands(1,:) = [0 fb(1)]; % all lower frequencies
		for b = 2:nbb
			fbands(b,:) = [fb(b-1) fb(b)];
		end
		fbands(nbb+1,:) = [fb(nbb) fs/2]; % up to Nyqvist frequency
	end
	nfbands = size(fbands,1);
end

% Sort out ROIs

nchans = size(V,1);
if ischar(rois) && strcmpi(rois,'allchans')    % all channels as a single ROI
	rois = {1:nchans};
	perchan = false;
	nrois = 1;
elseif ischar(rois) && strcmpi(rois,'perchan') % each channel as an ROI
	perchan = true;
	nrois = nchans;
else
	assert(isvector(rois) && iscell(rois),'ROI spec must be ''allchans'', ''perchan'', or a cell vector');
	perchan = false;
	nrois = length(rois);
end

% Calculate CPSD

S = var_to_cpsd(A,V,fres);
h = fres+1;

% Calculate log-determinants of CPSD

LDS = zeros(h,nrois);
if perchan
	for i = 1:nchans
		LDS(:,i) = log(squeeze(S(i,i,:))); % log-autopower!
	end
else
	for r = 1:nrois
		Sr = S(rois{r},rois{r},:);
		for k = 1:h
			LDS(k,r) = logdet(Sr(:,:,k));
		end
	end
end

if broadband
	erates = LDS;
else

	% Integrate across specified frequency bands

	erates = zeros(nfbands,nrois);
	for r = 1:nrois
		for b = 1:nfbands
			erates(b,r) =  bandlimit(LDS(:,r),[],fs,fbands(b,:));
		end
	end
%	fprintf('\nEntropy rates:\n');
%	disp(erates);
end
