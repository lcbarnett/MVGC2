function [A,C,K,V] = ss_normalise(A,C,K,V)

% Normalise ISS model by variance

G = ss_to_autocov(A,C,K,V,0); % covariance matrix

v = sqrt(diag(G));
u = 1./v;

C = u.*C;
K = K.*v';
V = u.*V.*u';
