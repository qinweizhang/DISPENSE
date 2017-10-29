function res = mtimes(a,b)
% res = mtimes(FT, x)
%
% reshape to 5D tensor of size (imsize(1),imsize(2),nc,s4,s5);
% ifft /fft
% reshape back to 1-unfolded tensor

F_for = @(x) fft2c(x);
F_adj = @(y) ifft2c(y);

S_adj = @(x) (sum(bsxfun(@times, conj(a.sens), x), 3));
S_for = @(x) bsxfun(@times, a.sens, x);

if(isempty(a.NUFFT_op))  
    %% nomal case, no mixed trajectory
    if size(b,1)==a.imsize(1)&& size(b,2)==a.imsize(2);
        res=b;
    elseif size(b,1)==a.imsize(1)*a.imsize(2)*a.ncoils || size(b,1)==a.imsize(1)*a.imsize(2)
        if a.adjoint
            res=reshape(b,[a.imsize(1),a.imsize(2),a.ncoils,numel(b)/(a.imsize(1)*a.imsize(2)*a.ncoils)]);
        else
            res=reshape(b,[a.imsize(1),a.imsize(2),1,numel(b)/(a.imsize(1)*a.imsize(2))]);
        end
        
    else
        error('unsupported matrix size!')
    end
    
    if a.adjoint
        res=S_adj(F_adj(res));
    else
        res=F_for(S_for(res));
    end
    
    
    if size(b,1)==a.imsize(1)*a.imsize(2)*a.ncoils || size(b,1)==a.imsize(1)*a.imsize(2)
        if a.adjoint
            res=reshape(res,[a.imsize(1)*a.imsize(2),numel(b)/(a.imsize(1)*a.imsize(2)*a.ncoils)]);
        else
            res=reshape(res,[a.imsize(1)*a.imsize(2)*a.ncoils,numel(b)/(a.imsize(1)*a.imsize(2))]);
        end
    end
    
    set_MCFop_adjoint(a,0); %temporary hack: reset adjoint value to 0 after every mtimes because we do not know how handle classes work :(
else
    %% SoSNav LRT recon: first column of kspace is 1D spiral data; second column of kspace is 2D cartesian TSE data
    if(size(b,1)==a.imsize(1)*a.imsize(2)*a.ncoils)  % b is in ksapce
        assert(a.adjoint,'b should be in the k-space!');
        
        k_data_4d = reshape(b,a.imsize(1)*a.imsize(2), a.ncoils, a.dimsize(1), a.dimsize(2));
        nav_k = k_data_4d(1:a.trj_length,:,:,1);
        tse_k = reshape(k_data_4d(:,:,:,2),a.imsize(1), a.imsize(2), a.ncoils, a.dimsize(1));
        
        % k-space --> image-space
        for s = 1:size(nav_k, 3)
            nav_i(:,:,1,s) = a.NUFFT_op' * col(nav_k(:,:,s));
        end
        tse_i = S_adj(F_adj(tse_k));
        
        res = cat(5, nav_i, tse_i);
        % reshape
        res=reshape(res,[a.imsize(1)*a.imsize(2),numel(res)/(a.imsize(1)*a.imsize(2))]);
        
    elseif(size(b,1)==a.imsize(1)*a.imsize(2)) % b is in image-space
        
        ima = reshape(b, a.imsize(1),a.imsize(2), 1, a.dimsize(1), a.dimsize(2));
        nav_i = ima(:,:,:,:,1);
        tse_i = ima(:,:,:,:,2);
        
        %image --> k-space
        for s = 1:size(nav_i, 4)
            nav_k(:,:,s) = reshape(a.NUFFT_op * nav_i(:,:,:,s), a.trj_length, a.ncoils);
        end
        tse_k=F_for(S_for(tse_i));
        
        %reshape: pad nav_k to cat with tse_k
        nav_k_pad = zeros(size(tse_k,1)*size(tse_k,2),size(tse_k,3),size(tse_k,4));
        nav_k_pad(1:size(nav_k,1),:,:) = nav_k; 
        nav_k_pad_2d = reshape(nav_k_pad, size(nav_k_pad,1)*size(nav_k_pad,2), size(nav_k_pad,3));
        tse_k_2d = reshape(tse_k, size(tse_k,1)*size(tse_k,2)*size(tse_k,3), size(tse_k,4));
        
        res = cat(2, nav_k_pad_2d, tse_k_2d);

        
    else
        error('Unsupported matrix size! kspace should be in size of [kx*ky*nc, dim1*dim2]; image should be in size of [kx*ky, dim1*dim2]');
    end
    
    set_MCFop_adjoint(a,0); %temporary hack: reset adjoint value to 0 after every mtimes because we do not know how handle classes work :(

    
end
