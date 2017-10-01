function C=conj_grad_G_2(G,C,B,Z,beta,du,Phi,F)
% new attempt
% minimizes the function
% argminG ||d - Fu G C Phi ||_2^2  + <Z,B-Z> - (beta/2) ||B -C ||_F ^2
disp('--- CG algorithm for C ---')

%params:
betals=0.6;
t0=1 ;
ls_alpha=1e-2;
maxlsiter=50;
maxiter=40;

iter=0;
s=0;
grad=calc_grad(G,F,du,Phi,B,C,beta); 

while(1)
    
    [f0,obj_l2,obj_inner_product,obj_F]= calc_objective(F,G,C,Phi,B,Z,du,beta);
    
    if iter==0
    disp(['iter: ',num2str(iter), '| obj: ',num2str(f0),'| obj_l2: ',num2str(obj_l2),'| obj_ip: ',num2str(obj_inner_product),'| obj_F: ',num2str(obj_F)])
    end
    
    % 1 calculate steepest direction
    gradk=calc_grad(G,F,du,Phi,B,C,beta);
    
    % 2 compute beta
    beta=(gradk(:).'*gradk(:))/(grad(:).'*grad(:)+eps); %Fletcher-Reeves;
    
    % 3 update conjugate direction
    sk=-gradk+ beta*s;
    
    % 4 perform a line search :
    t=t0;
    lsiter=0;
    f1=1e99; %only to start...
    while (f1 > f0 - ls_alpha*t*abs(gradk(:)'*sk(:)))^2 & (lsiter<maxlsiter)
        t = t * betals;
        [f1]  =   calc_objective(F,G,C+t*sk,Phi,B,Z,du,beta);
        lsiter=lsiter+1;
    end
    iter=iter+1;
    
    if lsiter>=maxlsiter | (norm(s(:)) < 1e-8) |iter==maxiter 
        disp('CG convergence reached.../ max line search/ max iters')
        break
    end
    
    % update the position
    C=C+t*sk;
    
    % print some iteration comments
    disp(['iter: ',num2str(iter),'| lsiter: ',num2str(lsiter), '| obj: ',num2str(f1),'| obj_l2: ',num2str(obj_l2),'| obj_ip: ',num2str(obj_inner_product),'| obj_F: ',num2str(obj_F)])

    %update parameters for next iteration;
    if lsiter > 2
		t0 = t0 * betals;
	end 
	if lsiter<1
		t0 = t0 / betals;
	end

    grad=gradk;
    s=sk;    
end
return


    function [obj,obj_l2,obj_inner_product,obj_F] = calc_objective(F,G,C,Phi,B,Z,du,beta)        
        obj_l2_inner= du - (F*(G*C*Phi)) ;
        obj_l2_inner=obj_l2_inner.*(abs(du)>0); % make this more efficientS
        obj_l2=obj_l2_inner(:)'*obj_l2_inner(:);
        
        obj_inner_product=trace(Z'*(B-Z));
        
        obj_F=norm((B-C),'fro')^2 ;
        
        obj=obj_l2+obj_inner_product+(beta/2).*obj_F;
    end

    function gradk=calc_grad(G,F,d,Phi,B,C,beta)
        gl2=grad_l2(G,F,d,Phi,C);
%         gip=grad_inner_product();
        gf=grad_F(B,C,beta);
        gradk=gl2+gf; %+gip+gf; %to do
    end


    function grad=grad_l2(G,F,d,Phi,C)
        grad=2*G'*(F'*((d+F*(G*C*Phi))*Phi')); % really dont know...
%         grad=-2*G'*(F'*((d+F*(G*C*Phi))*Phi')); % really dont know... (minus??)

    end

    function grad = grad_inner_product()
        % to do
    end

    function grad= grad_F(B,C,beta)
        grad= beta*(2*B+C);     %GUESS AT THIS POINT
    end

end
