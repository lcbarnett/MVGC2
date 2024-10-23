function V = var_predict(A,X)

% Return residuals covariance matrix of VAR prediction

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficient matrices not square');
A = reshape(A,n,n*p);   % concatenate VAR coefficient matrices

[n,m,N] = size(X);
assert(p < m,'too many lags or bad model order (p = %d, m = %d)',p,m);
M = N*(m-p); % effective number of observations

X = demean(X);          % remove temporal mean

% store lags

obs = p+1:m;
X0 = reshape(X(:,obs,:),n,M);
XL = zeros(n,p,M);
for k = 1:p
	XL(:,k,:) = reshape(X(:,obs-k,:),n,M);
end
XL = reshape(XL,p*n,M); % stack lagged observations

E = X0-A*XL;            % residuals
V = (E*E')/(M-1);       % residuals covariance matrix
