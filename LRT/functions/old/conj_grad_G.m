function Gk=conj_grad_G(G,C,Ak,Bk,Y,Z,alpha,Psi,du,Phi,P0_1)
% Gamma_k = Omega(F G C_k Phi)
% Gamma*_k(d) should go from measured data to estimate of G....
% Omega : row selection operator??
% Psi: Wavelet operator

% in paper; conjugate gradient, but this is a linear least squares??? 


id=(alpha/2).*(Psi'*(Ak+Y./alpha)); %size (id) = res^2 x Lg (size of G) 

% G_hat= 

% how to do the adjoint of Omega??

disp('how to find G_hat???')

sum(G_hat(:))

nom=G_hat+id;


PsiPsi=Psi*(Psi')*ones(size(id)); %% really at a loss here....
LambdaLambda=ones(size(id)); % no idea about this one...S

denom=LambdaLambda+(alpha/2).*PsiPsi; % TO DO...
denom=denom.^(-1);

disp('to do: denominator ')
Gk=denom.*nom;
end