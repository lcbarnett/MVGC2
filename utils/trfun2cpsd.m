function S = trfun2cpsd(H,V,fres,autospec)

if nargin < 4 || isempty(autospec), autospec = false; end

n = size(H,1);
h = fres+1;

S = NaN; % ensure a "bad" return value if anything goes wrong (see routine 'isbad')

[VL,cholp] = chol(V,'lower');
if cholp, return; end % show stopper

if autospec
	S = zeros(h,n);
	for k = 1:h
		HVLk = H(:,:,k)*VL;
		Sk = HVLk*HVLk';
		S(k,:) = diag(Sk);
	end
else
	S = zeros(n,n,h);
	for k = 1:h
		HVLk = H(:,:,k)*VL;
		S(:,:,k) = HVLk*HVLk';
	end
end
