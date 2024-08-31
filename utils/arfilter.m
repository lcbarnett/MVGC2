function Y = arfilter(A,X,allt)

% Y = arfilter(A,X,allt)
%
% Implements a vector AR filter with VAR coefficients A and input X.
%
% Implements:
%
%     Y(t) = X(t) +  sum_k A(k)*Y(t-k)
%
% if the flag allt is set, the above starts at t = 1, and the sum over k is for all
% indices for which the quantities are defined. If it is not set, the above is defined
% for t > p, and the sum is from k = 1 to k = p (so initial values in X are copied
% unchanged to Y).
%
% This function is a wrapper for the MEX routine 'arfilter_mex'.

global have_arfilter_mex;

assert(isreal(X) && isa(X,'double') && ismatrix(X), 'Input must be a real, double-precision matrix')

[n,m] = size(X);

if isempty(A)
	Y = X;
	return
end

assert(isreal(A) && isa(A,'double') && (ismatrix(A) || ndims(A) == 3), 'VAR coefficients must be a real, double-precision vector, matrix or 3D array')
avec = isvector(A); % vector of coefficients - apply filter to each row of X
if avec
	p = length(A);
else
	assert(size(A,1) == n && size(A,2) == n, 'VAR coefficient blocks must match input size');
	p = size(A,3);
end


if nargin < 3
	allt = true;
else
	assert(isscalar(allt) && islogical(allt), '''allt'' flag must be a logical scalar');

end

if ~allt && m <= p
	Y = X;
	return
end

if have_arfilter_mex
	Y = arfilter_mex(A,avec,p,X,allt);
	return
end

Y = X;
if allt
	if avec
		for t = 1:m
			for k = 1:p
				if k < t
					Y(:,t) = Y(:,t) + A(k)*Y(:,t-k);
				end
			end
		end
	else % ~avec
		for t = 1:m
			for k = 1:p
				if k < t
					Y(:,t) = Y(:,t) + A(:,:,k)*Y(:,t-k);
				end
			end
		end
	end
else % ~allt (can assume m > p)
	if avec
		for t = p+1:m
			for k = 1:p
				Y(:,t) = Y(:,t) + A(k)*Y(:,t-k);
			end
		end
	else % ~avec
		for t = p+1:m
			for k = 1:p
				Y(:,t) = Y(:,t) + A(:,:,k)*Y(:,t-k);
			end
		end
	end
end
