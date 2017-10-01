function Gk=conj_grad_G_3(G,C,A,Y,alpha,Psi,d,Phi,F)
% TO DO: THINK VERY WELL ABOUT MASK (OMEGA) 
tol=10;
maxiter=10;

L=Lambda(F,C,Phi,abs(d)>0);
a2PP=(alpha/2)*Psi'*Psi;

b=(L'*d + (alpha/2)*Psi'*(A+Y./alpha));
X= @(G) (L'*(L*G)) +a2PP*G;

%% TEMP: USE SPOT OPERATORS
%{
1+1
bspot=reshape(b,[numel(b) 1]);
Xspot= @(G,LL) (L'*(L*G)) +a2PP*reshape(G,[numel(G)/LL,LL]); % should be some sort of blockdiagonal operator of Lrank ops
spotoperator=opFunction(6553600,6553600,Xspot,1,1)

solution=bicg(spotoperator,bspot)
% Xa= @(d) (L*(L'*d)) +a2PP'*d; %not needed
%}
%%
x=G;
R=b-X(x);
P=R;


iter=0;resid=1e99;
while iter<maxiter && resid>tol
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
Gk=xk;

end
