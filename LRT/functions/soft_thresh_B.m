function [Bk,mu]=soft_thresh_B(C,Z,mu_o,beta,params)
mu=mu_o;

if params.iter>1 || params.automu==0
    mub=mu/beta;
    CZb=C-(Z./beta); %precalculate to speed up
    thr=abs(CZb)>=mub; %threshold
    
else        %first iter and automu - 
    fprintf('Automatically selecting mu...\n')
    
    nthr=0;
    iter_automu=0;
    while (nthr<0.3)||(nthr>0.9) 
        iter_automu=iter_automu+1; 
        
        if nthr<0.3;
            mu=mu.*0.9;
        elseif nthr>0.9
            mu=mu.*1.1;
        end
        if iter_automu>100
        mu =mu_0; 
        fprintf('automu failed \n')
        end
        
        mub=mu/beta;
        CZb=C-(Z./beta); %precalculate to speed up
        thr=abs(CZb)>=mub; %threshold
        nthr=sum(thr(:))./numel(thr);
        if iter_automu>100
            break;
        end
    end
            fprintf('mu= %d \n',mu)
end

fprintf('soft-thresholding B: ')
fprintf('mu/beta: %1.2e  |',(mub))
fprintf('thresholding %1.2f  percent \n',(100*(sum(~thr(:))./numel(thr))))


Bk=thr.*(abs(CZb)-mub).*((CZb)./(abs(CZb)+eps));

title_text = sprintf('s.t. C: %d data points thresholded (%1.2f%%).',sum(~thr(:)), (100*(sum(~thr(:))./numel(thr))));
figure(998);subplot(222); spy(thr); title(title_text); 
end