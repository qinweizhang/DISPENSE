function res = mtimes(a,b)
% res = mtimes(FPS, x)


F_for = @(x) fft2d(x);      %2D forward FFT
F_adj = @(x) ifft2d(x);     %2D backward FFT

S_adj = @(x) (sum(bsxfun(@times, conj(a.sens), x), 3));     %a.sens: in [ny, nz, nc]; x in [ny, nz, nc]; output: [ny nz]
S_for = @(x) bsxfun(@times, a.sens, x);                     %a.sens: in [ny, nz, nc]; x in [ny, nz]; output: [ny nz, nc]

P_adj = @(x) (sum(bsxfun(@times, conj(a.pe), x), 4));       %a.pe in [ny, nz, 1, nshot]; x in [ny, nz, nc, nshot]; output: [ny, nz, nc]
P_for = @(x) bsxfun(@times, a.pe, x);                       %a.pe in [ny, nz, 1, nshot]; x in [ny, nz, nc]; output: [ny, nz, nc, nshot]



if size(b,1)==prod(a.image_dim)*a.numCoils*a.numShots || numel(b)==prod(a.image_dim)
    if a.adjoint
        res=reshape(b,[a.image_dim(1),a.image_dim(2), a.numCoils,a.numShots]);
    else
        res=reshape(b,[a.image_dim(1),a.image_dim(2)]);
    end
    
else
    error('unsupported matrix size!')
end



if a.adjoint
    res=S_adj(P_adj(F_adj(res)));
else
    res=F_for(P_for(S_for(res)));
    res=res.*a.mask;
    res = col(res);
end




