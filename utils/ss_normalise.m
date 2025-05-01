function [A,C,K,V] = ss_normalise(A,C,K,V)

% Transform ISS model so that process covariance matrix = I.
%
% NOTE: sub-process covariance matrices will also be identity.

% Construct the transformation

G = ss_to_autocov(A,C,K,V,0); % process covariance matrix
U = chol(G,'lower');          % the inverse transformation
T = inv(U);                   % the transformation

% Transform the ISS parameters

C = T*C;
K = K*U;
L = T*chol(V,'lower');
V = L*L';
