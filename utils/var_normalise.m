function [A,V] = var_normalise(A,V)

% Transform VAR model so that process covariance matrix = I.
%
% NOTE: sub-process covariance matrices will also be identity.

% Construct the transformation

G = var_to_autocov(A,V,0); % process covariance matrix
U = chol(G,'lower');       % the inverse transformation
T = inv(U);                % the transformation

% Transform the VAR parameters

for k = 1:size(A,3);
	A(:,:,k) = T*A(:,:,k)*U;
end
L = T*chol(V,'lower');
V = L*L';
