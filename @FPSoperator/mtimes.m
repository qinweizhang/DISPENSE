function res = mtimes(a,b)
% res = mtimes(FPS, x)


F_for = @(x) fft3d(x);      %3D forward FFT
F_adj = @(x) ifft3d(x);     %3D backward FFT

S_adj = @(x) (sum(bsxfun(@times, conj(a.sens), x), 4));     %a.sens: in [nx, ny, nz, nc]; x in [nx, ny, nz, nc]; output: [nx ny nz]
S_for = @(x) bsxfun(@times, a.sens, x);                     %a.sens: in [nx, ny, nz, nc]; x in [nx, ny, nz]; output: [nx ny nz, nc]

P_adj = @(x) (sum(bsxfun(@times, conj(a.pe), x), 5));       %a.pe in [nx, ny, nz, 1, nshot]; x in [nx, ny, nz, nc, nshot]; output: [nx, ny, nz, nc]
P_for = @(x) bsxfun(@times, a.pe, x);                       %a.pe in [nx, ny, nz, 1, nshot]; x in [nx, ny, nz, nc]; output: [nx, ny, nz, nc, nshot]



if size(b,1)==prod(a.image_dim)*a.numCoils*a.numShots || size(b,1)==prod(a.image_dim)
    if a.adjoint
        res=reshape(b,[a.image_dim(1),a.image_dim(2),a.image_dim(3), a.numCoils,a.numShots]);
    else
        res=reshape(b,[a.image_dim(1),a.image_dim(2),a.image_dim(3)]);
    end
    
else
    error('unsupported matrix size!')
end



if a.adjoint
    res=S_adj(P_adj(F_adj(res)));
else
    res=F_for(P_for(S_for(res)));
    res=res.*a.mask;
end


res = col(res);


