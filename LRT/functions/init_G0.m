function [Phi,G0,C0,A,B,Y,Z]=init_G0(P,Psi,nav_estimate_1,nav_estimate_2,L1)
%initializes all matrices.

%input: 
% P: 1-unfolded zero-filled recon
% Phi: kronecker product of subspace estimates
% L1: estimated rank of G
disp('initializing matrices...')
Phi=kron(nav_estimate_2,nav_estimate_1).';      %from subspaces Phi= kron(G^4,G^3)^T

X=P*Phi'; 
% [U,S,V] = svds(X,L1);
[U,S,V] = svd(X, 'econ'); 

G0=U(:,1:L1);                           % first Lg vectors from left-dominant  svd of P10 Psi^H

C0=G0'*P*Phi';                          % G0^H P10 Psi^H 

A=zeros(size(Psi*G0));
Y=zeros(size(A));

B = zeros(size(C0));
Z = zeros(size(C0));
end
