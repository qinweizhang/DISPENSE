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

