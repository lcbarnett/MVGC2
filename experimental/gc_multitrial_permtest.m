function [F1,F2,pval,varmo1,varmo2,DFp] = gc_multitrial_permtest(X1,X2,x,y,nperms,varmomax,varmosel,varmoreg,verb,fignum,nhist)

% A permutation test for multitrial data, for testing whether a GC value is different between conditions.
%
% NOTE: The null hypothesis here isn't, strictly speaking, that the GCs are different, but rather that
% the data in the two conditions is drawn from a different distribution. In practice, however, it seems
% unlikely that if the distributions are different that the GCs *won't* be different! Thus a rejection
% of the null in this case still tells us something useful.
%
% This routine only handles VAR estimation; in principle, we can do the same for state-space estimation.
%
% INPUTS
%
% X1         - multitrial time series data for condition 1
% X2         - multitrial time series data for condition 2
% x          - vector of indices of target variable
% y          - vector of indices of source variable (all remaining variables are conditioned on)
% nperms     - number of samples for the empirical permutation null distribution
% varmosel   - VAR model order selection criterion ('AIC', 'BIC', 'HQC', 'LRT')
% varmomax   - maximum VAR model order for model order estimation
% varmoreg   - regression method for model order estimation ('OLS' or 'LWR')
% verb       - verbosity flag
% fignum     - figure number for histogram display (set to 0 for no display)
% nhist      - number of bins for the histogram display
%
% OUTPUTS
%
% F1         - GC estimate for condition 1
% F2         - GC estimate for condition 2
% pval       - p-value for 2-tailed test of F2-F1 against the permutation null distribution
% varmo1     - VAR model order estimate for condition 1
% varmo2     - VAR model order estimate for condition 2
% DFp        - the empirical permutation null distribution

if nargin <  8 || isempty(varmoreg),  varmoreg  = 'OLS'; end
if nargin <  9 || isempty(verb),      verb      = true;  end
if nargin < 10 || isempty(fignum),    fignum    = 0;     end
if nargin < 11 || isempty(nhist),     nhist     = 40;    end

[nvars1,nobs1,ntrials1] = size(X1);
[nvars2,nobs2,ntrials2] = size(X2);

assert(nvars2 == nvars1 && nobs2 == nobs1,'Data series don''t match!');

% Model condition 1 data and calculate sample GC

[moaic,mobic,mohqc,molrt] = tsdata_to_varmo(X1,varmomax,varmoreg,[],[],[],verb);
varmo1 = selvarmo(varmosel,moaic,mobic,mohqc,molrt);
[A1,V1] = tsdata_to_var(X1,varmo1,'OLS');
F1 = var_to_mvgc(A1,V1,x,y);

% Model condition 1 data and calculate sample GC

[moaic,mobic,mohqc,molrt] = tsdata_to_varmo(X2,varmomax,varmoreg,[],[],[],verb);
varmo2 = selvarmo(varmosel,moaic,mobic,mohqc,molrt);
[A2,V2] = tsdata_to_var(X2,varmo2,'OLS');
F2 = var_to_mvgc(A2,V2,x,y);

DF = F2-F1; % difference in GCs

if verb	% report
	RDF = 2*DF/(F1+F2); % relative magnitude of GC difference
	fprintf('F1 (estimate) = % 8.6f : varmo(%s) = %2d\n',F1,varmosel,varmo1);
	fprintf('F2 (estimate) = % 8.6f : varmo(%s) = %2d\n',F2,varmosel,varmo2);
	fprintf('DF (estimate) = % 8.6f (%8.6f)\n',DF,RDF);
end

% Pool trials

Xp = cat(3,X1,X2);
ntrialsp = ntrials1+ntrials2;

% Permutation test - we use the original model orders

if verb
	fprintf('\nPermutation test ')
	k10 = round(nperms/10);
	tstart = tic;
end

DFp = zeros(nperms,1);
for k = 1:nperms

	sidx = randperm(ntrialsp); % shuffled indices for pooled trials

	% Split shuffled pooled sample, model, and calculate GCs

	X1p = Xp(:,:,sidx(1:ntrials1));
	[A1,V1] = tsdata_to_var(X1p,varmo1,'OLS');
	F1p = var_to_mvgc(A1,V1,x,y);

	X2p = Xp(:,:,sidx(ntrials1+1:end));
	[A2,V2] = tsdata_to_var(X2p,varmo2,'OLS');
	F2p = var_to_mvgc(A2,V2,x,y);

	% The observed difference in permutation sample GCs

	DFp(k) = F2p-F1p;

	if verb && rem(k,k10) == 0, fprintf('.'); end % progress indicator

end

if verb
	tend = toc(tstart);
	dur = seconds(tend);
	dur.Format = 'hh:mm:ss.SS';
	durstr = regexprep(char(dur),'^00:(00:)?','');
	fprintf(' %s (hh:mm:ss.SS)\n',durstr);
end

% 2-tailed p-value

pval = sum(abs(DFp) >= abs(DF))/nperms;
if verb, fprintf('\np-value = %f\n\n',pval); end

if fignum > 0

	% Display histogram of permutation null distribution of differences

	figure(fignum); clf
	DFpm = median(DFp);
	col = [0.2 0.4 0.8];
	histogram(DFp,nhist,'FaceColor',col);
	xline(DF,  'r','LineWidth',2,'Label','DF');
	xline(DFpm,'g','LineWidth',2,'Label','DFp median');
	title('GC Permutation test');
	xlabel('DFp (permutation null)');
	ylabel('Frequency');
	grid on

end

end

function varmo = selvarmo(varmosel,moaic,mobic,mohqc,molrt)

	switch upper(varmosel)
		case 'AIC', varmo = moaic;
		case 'BIC', varmo = mobic;
		case 'HQC', varmo = mohqc;
		case 'LRT', varmo = molrt;
		otherwise, error ('Unknown model order selection criterion');
	end

end
