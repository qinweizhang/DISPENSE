function res = mtimes(a,b)
% res = mtimes(FPS, x)

% S = SUM(S(l));  S(l) = interpo_l .* omiga_l * NUFFT * I


if a.adjoint %kspa --> ima
    %%
    assert(length(b)/a.ch_nr==size(a.interpolator ,2),'interpolator size doesnot match with data');

    sig_ch_by_ch = reshape(b, size(a.interpolator ,2), 1,  a.ch_nr);
    
    S_tseg = bsxfun(@times, sig_ch_by_ch, conj(a.interpolator)');   %S_seg now with weighting from interpolator in size of [samples, segs, ch_nr]
    
    I_tseg = zeros([a.image_dim, (a.tSeg.nr_segments+1)]);
    for k = 0:a.tSeg.nr_segments
        I_tseg(:,:,:,k+1) = (a.nuFTop' * col(S_tseg(:,k+1,:))) .* exp(1i * a.b0_maps_mHz * a.tSeg.tau * k);
    end
    res = sum(I_tseg, 4);
    
    
    
    
else %ima --> kspa
    %%
    S_tseg_no_weighting = zeros(length(a.samples_idx),(a.tSeg.nr_segments+1), a.ch_nr);  %S_seg now without weighting from interpolator in size of [samples, segs, ch_nr]
    for k = 0:a.tSeg.nr_segments
        S_tseg_no_weighting(:,k+1,:) =  reshape( (a.nuFTop *  (b .* exp(-1i * a.b0_maps_mHz * a.tSeg.tau * k))), length(a.samples_idx), a.ch_nr);       
    end
    
    
    %apply interpolater
    S_tseg = bsxfun(@times, S_tseg_no_weighting, a.interpolator');
    
    res = col(sum(S_tseg, 2));
    
    
end









