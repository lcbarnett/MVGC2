function [A,V] = var_studentise(A,V)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[nn1,nn2] = size(V);
assert(nn1 == nn2,'residuals covariance matrix not square');
assert(nn1 == n  ,'residuals covariance matrix doesn''t match VAR coefficients matrix');

% Normalise VAR model by process variance

G = var_to_autocov(A,V,0)

% TODO - normalisation should make G = I !!
