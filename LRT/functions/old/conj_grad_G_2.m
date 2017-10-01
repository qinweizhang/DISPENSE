function G=conj_grad_G_2(G,C,A,Y,alpha,Psi,du,Phi,F)
% new attempt
% minimizes the function
% argminG ||d - Fu G C Phi ||_2^2  + <Y,A-Psi G> + (alpha/2) ||A - Psi G ||_F ^2
% argminG ||d - Fu G C Phi ||_2^2  + <Y,A> -<Y,Psi G> + (alpha/2) ||A - Psi G ||_F ^2
% TO DO: CHECK GRADIENTS
% UNCOMMENT GRADIENT OF INNER PRODUCT 
% TUNE PARAMS
% MORE TESTING
% OPTIMIZE SPEED WITH PRECALCULATING FFTS AND ADDING SMARTER MASKING IN OBJ
% FUNTCTION 
disp('--- CG algorithm for G ---')

%params:
betals=0.6;
t0=1 ;
ls_alpha=1e-2;
maxlsiter=50;
maxiter=40;


iter=0;
s=0;
grad=calc_grad(G,C,Phi,du,F,Y,Psi,alpha); 

while(1)
    [PsiG, Psis, FGCPhi, FsCPhi]=preobjective(F,G,C,Phi,Psi,grad);
    [f0,obj_l2,obj_inner_product,obj_F]= calc_objective(A,Y,PsiG, Psis, FGCPhi, FsCPhi,0,alpha,du);
    
    if iter==0
    disp(['iter: ',num2str(iter), '| obj:',num2str(f0),'| obj_l2: ',num2str(obj_l2),'| obj_ip: ',num2str(obj_inner_product),'| obj_F: ',num2str(obj_F)])
    end
    
    % 1 calculate steepest direction
    gradk=calc_grad(G,C,Phi,du,F,Y,Psi,alpha);
    
    % 2 compute beta
    beta=(gradk(:).'*gradk(:))/(grad(:).'*grad(:)+eps); %Fletcher-Reeves;
    
    % 3 update conjugate direction
    sk=-gradk+ beta*s;
    
    % 3b : calculate preobjective
    [PsiG, Psis, FGCPhi, FsCPhi]=preobjective(F,G,C,Phi,Psi,sk);
    
    % 4 perform a line search :
    t=t0;
    lsiter=0;
    f1=1e99; %only to start...
    while (f1 > f0 - ls_alpha*t*abs(gradk(:)'*sk(:)))^2 & (lsiter<maxlsiter)
        t = t * betals;
        [f1,obj_l2,obj_inner_product,obj_F]  =   calc_objective(A,Y,PsiG, Psis, FGCPhi, FsCPhi,t,alpha,du);
        lsiter=lsiter+1;
    end
    iter=iter+1;
    
    if lsiter>=maxlsiter | (norm(s(:)) < 1e-8) | iter==maxiter
        disp('CG convergence reached.../ max line search/ max iters')
        break
    end
    
    % update the position
    G=G+t*sk;
    
    % print some iteration comments
    disp(['iter: ',num2str(iter),'| lsiter: ',num2str(lsiter), '| obj: ',num2str(f1),'| obj_l2: ',num2str(obj_l2),'| obj_ip: ',num2str(obj_inner_product),'| obj_F: ',num2str(obj_F)])

    %update parameters for next iteration;
    if lsiter > 2
		t0 = t0 * betals;
	end 
	
	if lsiter<1
		t0 = t0 / betals;
	end

    s=sk;
    grad=gradk;

    
end
return

    function [PsiG, Psis, FGCPhi, FsCPhi]=preobjective(F,G,C,Phi,Psi,sk); %so we only have to calculate this once per outer iteration
        PsiG=Psi*G;
        Psis=Psi*sk;
        FGCPhi=F*(G*C*Phi);
        FsCPhi=F*(sk*C*Phi);
    end

    function [obj,obj_l2,obj_inner_product,obj_F] = calc_objective(A,Y,PsiG, Psis, FGCPhi, FsCPhi,t,alpha,du)
        %objective= ||d - Fu G C Phi ||_2^2  + <Y,A-Psi G> + (alpha/2) ||A - Psi G ||_F ^2      
        obj_l2_inner= du - (FGCPhi+t*FsCPhi) ;
        obj_l2_inner=obj_l2_inner.*(abs(du)>0); % make this more efficient
        obj_l2=obj_l2_inner(:)'*obj_l2_inner(:);
        
        obj_inner_product=trace(Y'*(A-(PsiG+t*Psis))); %additive
        
        obj_F=norm((A-(PsiG+t*Psis)),'fro')^2 ;        %sub additive
        
        obj=obj_l2+obj_inner_product+(alpha/2).*obj_F;
    end

    function gradk=calc_grad(G,C,Phi,du,F,Y,Psi,alpha)
        gl2=grad_l2(G,C,Phi,du,F);
        gip=grad_inner_product(Y,Psi);
        gf=grad_F(alpha,A,G,Psi);
        %         gradk=gl2+gip+gf; %temp: remove gip
        gradk=gl2+gf;

    end

    function grad=grad_l2(G,C,Phi,du,F)
        grad = 2* (F'*(du * Phi'*C' + F*(G*C*Phi*Phi'*C'))); %new attempt...
    end

    function grad = grad_inner_product(Y,Psi)
        % <y, Psi * G>
        grad= -2*(Y'*Psi); % not sure at all..
    end

    function grad= grad_F(alpha,A,G,Psi)
        grad= alpha*(2*Psi'*A+G); 
    end

end
