function Ck=conj_grad_C_3(G,C,B,Z,beta,d,Phi,F)
% TO DO: THINK VERY WELL ABOUT MASK (OMEGA) 
tol=1e-18;
maxiter=30;

L=Lambda2(F,G,Phi,abs(d)>0);

b=((L'*d) + (beta/2)*(B+Z./beta));
X= @(G) (L'*(L*G)) +(beta/2)*G;

%%
x=C;
R=b-X(x);
P=R;

iter=0;resid=1e99;
while iter<maxiter && resid > tol
    iter=iter+1;
    alpha= trace(R'*R)/trace(P'*X(P));
    xk=x+alpha.*P;
    Rk=R-alpha.*(X(P));
    
    resid=abs(sum(Rk(:)));
    disp(['iteration: ',num2str(iter),'|residual: ',num2str(resid)])
    
    beta=trace(Rk'*Rk)/trace(R'*R);
    Pk=Rk+beta.*P;
    R=Rk; P=Pk; x=xk;
end
Ck=xk;





end
