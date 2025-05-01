function [KR,VR,rep] = vardarea(VARA,V,r)

% Solve DARE for SS model derived from reduced VAR this is a simplified
% version of vardare where the sources comprise the whole system.

[n,~,p] = size(VARA);
nr  = length(r);
pn  = p*n;
pn1 = pn-n;

A = [reshape(VARA,n,pn); eye(pn1) zeros(pn1,n)];
C = reshape(VARA(r,:,:),nr,pn);
Q = [V zeros(n,pn1); zeros(pn1,pn)];
S = [V(:,r); zeros(pn1,nr)];
R = V(r,r);
[KR,VR,rep] = mdare(A,C,Q,R,S);
