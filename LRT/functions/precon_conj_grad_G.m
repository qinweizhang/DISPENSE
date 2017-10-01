function Gk=precon_conj_grad_G(G,C,A,Y,alpha,Psi,d,Phi,F,params)
tic; 
samplingmask=abs(d)>0; 
tol=params.G.tol;
maxiter=params.G.maxiter;

L=Lambda(F,C,Phi,samplingmask);
% a2PP=(alpha/2)*pinv(Psi)*Psi;
a2PP=(alpha/2)*Psi'*Psi;

%input data
% b=(L'*d + (alpha/2)*pinv(Psi)*(A+Y./alpha));
b=(L'*d + (alpha/2)*Psi'*(A+Y./alpha));
b=b(:);


% try to reshape operator so we have proper matrix, vector calculations
X= @(G) (L'*((L*G))) +a2PP*G;
Res= @(x) reshape(x,[numel(G),1]);
ResA= @(x) reshape(x,size(G));
Aop = @(G) Res(X(ResA(G))); 

if params.G.precon
    fprintf('precon CG for G_k: ')

    %build preconditioner
    CPPC=(C*Phi*Phi'*C');
    w=diag(CPPC)+alpha/2;
    w=permute(w,[2 1]);
    w=1./w;
    
    mfun=@(x) Res(bsxfun(@times, ResA(x),w)); %test
else
    fprintf('CG for G_k: ')
    mfun=[];
end

[x,flag,relres,iter,resvec]=bicgstab(Aop,b,tol,maxiter,mfun,[],G(:)); %add initial guess
% [x,flag,relres,iter,resvec]=cgs(Aop,b,tol,maxiter,mfun,[],G(:)); %add initial guess
% [x,flag,relres,iter,resvec]=pcg(Aop,b,tol,maxiter,mfun,[],G(:)); %add initial guess

Gk=ResA(x);

figure(998);subplot(223);
plot(log10(resvec./norm(b(:))),'r*-'); 
xlabel('iterations'); ylabel('10log of relative residual')
t=toc; 
fprintf('t: %i seconds',t)
fprintf('| relres %d | iters: %i | \n',relres,iter)

end
