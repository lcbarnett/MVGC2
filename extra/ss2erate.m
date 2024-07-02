function erates = ss2erate(A,C,K,V,rois,normalise)

% Calculate time-domain entropy rates for SS model for specified ROIs.
%
% A, C, K, are the ISS parameter matrices, V the residuals covariance matrix
%
% rois is a cell vector, where each cell is a vector of channel
% numbers for the channels in an ROI. It may also take two special
% values: 'allchans' treats the entire system as a single roi;
% 'perchan' returns entropy rates for each individual channel.
%
% If the normalise flag is set to true, then the ISS model is normalised
% for unit variance. Error rates are then equivalent to minus the mutual
% information between the process and its own past and as such are scale-
% invariant (and non-positive).
%
% Results are returned in the vector erates.

if nargin < 5 || isempty(rois) || (ischar(rois) && strcmpi(rois,'allchans')) % all channels as a single ROI (no DARE required)
	erates = logdet(V);
	return
end

if nargin < 6 || isempty(normalise)
	normalise = false;
end

if normalise
	[A,C,K,V] = ss_normalise(A,C,K,V);
end

% Calculate log-determinants of reduced residuals covariance matrices

[L,cholp] = chol(V,'lower');
assert(cholp == 0,'Residuals covariance matrix not positive-definite!');
KL = K*L;
KLK = KL*KL';

if ischar(rois) && strcmpi(rois,'perchan') % each channel as an ROI

	nchans = size(C,1);
	erates = zeros(nchans,1);
	for i = 1:nchans
		[~,VR,rep] = mdare(A,C(i,:),KLK,V(i,i),K*V(:,i)); % reduced model residuals covariance
		if sserror(rep,i), continue; end                  % check DARE report, bail out on error
		erates(i) = log(VR);                              % scalar, so no determinant required!
	end

else

	assert(isvector(rois) && iscell(rois),'ROI spec must be ''allchans'', ''perchan'', or a cell vector');
	nrois = length(rois);
	erates = zeros(nrois,1);
	for r = 1:nrois
		i = rois{r};
		[~,VR,rep] = mdare(A,C(i,:),KLK,V(i,i),K*V(:,i)); % reduced model residuals covariance
		if sserror(rep,r), continue; end                  % check DARE report, bail out on error
		erates(r) = logdet(VR);
	end

end

%fprintf('\nEntropy rates:\n');
%disp(erates);
