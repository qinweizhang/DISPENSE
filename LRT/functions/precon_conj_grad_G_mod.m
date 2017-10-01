function Gk=precon_conj_grad_G_mod(G,C,A,Y,alpha,Psi,d,Phi,F)
tic; 
fprintf('CG for G_k: ')

tol=1e-13;
maxiter=50; %temp

% L=Lambda(F,C,Phi,abs(d)>0);
a2PP=(alpha/2)*Psi'*Psi;

%input data
b=(((F'*d)*(Phi'*C')) + (alpha/2)*Psi'*(A+Y./alpha));
b=b(:);

%build preconditioner
CPPC=(C*Phi*Phi'*C');
w=diag(CPPC)+alpha/2; 
w=permute(w,[2 1]);
w=1./w;


% try to reshape operator so we have proper matrix, vector calculations
AA=(F*(G*C));
mask=double(abs(d)>0);


% Q=Phi(:,1)*Phi(:,1).';
X1= @(G) (F*(G*C));
X1b = @(G) loopfunc(X1(G),Phi,mask,C);
X2= @(G) F'*(X1b(G));
X= @(G) X2(G) +a2PP*G;
% X= @(G) (L'*(L*G)) +a2PP*G;
Res= @(x) reshape(x,[numel(G),1]);
ResA= @(x) reshape(x,size(G));
Aop = @(G) Res(X(ResA(G))); 
mfun=@(x) Res(bsxfun(@times, ResA(x),w)); %test

[x,flag,relres,iter,resvec]=cgs(Aop,b,tol,maxiter,[],[],G(:)); %add initial guess


Gk=ResA(x);
figure(28); plot(log10(resvec./norm(b(:))),'r*-'); 
xlabel('iterations'); ylabel('10log of relative residual')
t=toc; 
fprintf('Time taken: %d seconds',t)
fprintf('| relres %d | iters: %d | \n',relres,iter)

end

