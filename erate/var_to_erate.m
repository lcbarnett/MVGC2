function erates = var_to_erate(A,V,rois)

% Calculate time-domain entropy rates for SS model for specified ROIs.
%
% A is the VAR coefficients array, V the residuals covariance matrix
%
% rois is a cell vector, where each cell is a vector of channel
% numbers for the channels in an ROI. It may also take two special
% values: 'global' treats the entire system as a single roi;
% 'perchan' returns entropy rates for each individual channel.
%
% Results are returned in the vector erates.

if ischar(rois) && strcmpi(rois,'global') % all channels as a single ROI (no DARE required)
	erates = logdet(V);
	return
end

nchans = size(A,1);

if ischar(rois) && strcmpi(rois,'perchan') % each channel as an ROI

	erates = zeros(nchans,1);
	for i = 1:nchans
		[~,VR,rep] = vardarea(A,V,i);    % reduced model residuals covariance
		if sserror(rep,i), continue; end % check DARE report, bail out on error
		erates(i) = log(VR);             % scalar, so no determinant required!
	end

else

	assert(isvector(rois) && iscell(rois),'ROI spec must be ''global'', ''perchan'', or a cell vector');
	nrois = length(rois);
	erates = zeros(nrois,1);
	for r = 1:nrois
		i = rois{r};
		[~,VR,rep] = vardarea(A,V,i);    % reduced model residuals covariance
		if sserror(rep,r), continue; end % check DARE report, bail out on error
		erates(r) = logdet(VR);
	end

end

%fprintf('\nEntropy rates:\n');
%disp(erates);
