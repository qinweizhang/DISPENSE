%
%
% INPUT
%
% k_spa:              cpx phase inhoerent k space data in [ky kz n_coil n_shots], navigator is part of it.
% trj:                trajectory for k_spa. TODO
% sensemaps:          cpx sense maps in [ky kz n_coil]
% pars:               reconstruction parameters structure
%
% OUTPUT
% sol:                solution in [ky kz]
%
% (c) q.zhang 2017 Amsterdam

function sol = POCS_SENESE(kspa, trj, sense_map, pars)
[ky kz n_coil n_shots] = size(kspa);
x = zeros(ky, kz, n_shots); %initialize
x_esti = zeros(ky, kz);
it = 0;

W = zeros(ky, kz);
ky_range = floor(ky/2-pars.POCS.kernelsize(1)/2):floor(ky/2+pars.POCS.kernelsize(1)/2)-1;
kz_range = floor(kz/2-pars.POCS.kernelsize(2)/2):floor(kz/2+pars.POCS.kernelsize(2)/2)-1;
W(ky_range, kz_range) = 1;


while(it<pars.POCS.nit || residual_error>pars.POCS.tol)
    
    for sh=1:n_shots
        px_allshots(:,:,sh) = POCS_projection(x(:,:,sh), sense_map, kspa(:,:,:,sh), trj);
    end
    
    %shot-average
    fpx_allshots = fft2d(px_allshots);
    fpx_allshots_lf = bsxfun(@times, fpx_allshots, W);
    px_allshots_lr = ifft2d(fpx_allshots_lf);
    
    px_aver = mean(px_allshots.* (conj(px_allshots_lr) ./ abs(px_allshots_lr)),3); %remove lowresolution phase
    
    x_esti_old = x_esti;
    x_esti = x_esti_old + pars.POCS.lamda.*(px_aver - x_esti_old);      %update extimation
    x = bsxfun(@times, x_esti, (px_allshots_lr./abs(px_allshots_lr)));  %update initial point for next step (recover phase)
    
    residual_error = trace((x_esti-x_esti_old)'*(x_esti-x_esti_old)) / trace(x_esti_old'*x_esti_old);  
    
    it = it +1;
    
    disp(['nit: ', num2str(it)]); 
    disp(['residual_error: ', num2str(residual_error)]); 
    
    figure(11); subplot(121); imshow(abs(x_esti),[]); subplot(122); imshow(angle(x_esti),[-pi pi]);
    
    
    
end

sol =  x_esti;

end