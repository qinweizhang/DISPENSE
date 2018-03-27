function res = mtimes(a,b)
% res = mtimes(A_opt, x)
% measured_kspa in size of [ny nz  nc  ns]
% reconed_kspa in size of  [ny nz  ns]


F_for = @(x) fft2d(x);      %2D forward FFT
F_adj = @(x) ifft2d(x);     %2D backward FFTc

S_adj = @(x) (sum(bsxfun(@times, conj(a.sens), x), 3));     %a.sens: in [ny, nz, nc]; x in [ny, nz, nc, ns]; output: [ny, nz, ns]
S_for = @(x) bsxfun(@times, a.sens, x);                     %a.sens: in [ny, nz, nc]; x in [ny, nz, ns]; output: [ny nz, nc, ns]


if numel(b)==prod(a.kspa_dim)*a.numCoils*a.numShots || numel(b)==(prod(a.image_dim) * a.numShots )
    if a.adjoint %measured_kspa --> reconed_kspa
        res=reshape(b,[a.kspa_dim(1),a.kspa_dim(2), a.numCoils,a.numShots]);
    else %reconed_kspa --> measured_kspa
        res=reshape(b,[a.image_dim(1),a.image_dim(2),1, a.numShots]);
    end
    
else
    error('unsupported matrix size!')
end



if a.adjoint  %measured_kspa --> reconed_kspa
    res=squeeze(F_for(S_adj(F_adj(res))));
else  %reconed_kspa --> measured_kspa
    res=F_for(S_for(F_adj(res)));
    res=bsxfun(@times, res, permute(a.mask, [1 2 4 3]));
    res = col(res);
end




