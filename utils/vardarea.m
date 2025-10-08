function [K,V,rep] = vardarea(VARA,VARV,r)

% Solve DARE for SS model derived from reduced VAR this is a simplified
% version of vardare where the sources comprise the whole system.

[n,~,p] = size(VARA);
nr  = length(r);
pn  = p*n;
pn1 = pn-n;

A = [reshape(VARA,n,pn); eye(pn1) zeros(pn1,n)];
C = reshape(VARA(r,:,:),nr,pn);
Q = [VARV zeros(n,pn1); zeros(pn1,pn)];
S = [VARV(:,r); zeros(pn1,nr)];
R = VARV(r,r);
[K,V,rep,L,P] = mdare(A,C,Q,R,S);
