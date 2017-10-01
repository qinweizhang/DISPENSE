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
it = 1;
residual_error = inf;

W = zeros(ky, kz);
ky_range = floor(ky/2-pars.POCS.Wsize(1)/2):floor(ky/2+pars.POCS.Wsize(1)/2)-1;
kz_range = floor(kz/2-pars.POCS.Wsize(2)/2):floor(kz/2+pars.POCS.Wsize(2)/2)-1;
W(ky_range, kz_range) = 1;


while(it<pars.POCS.nit && residual_error>pars.POCS.tol)
    
    for sh=1:n_shots
        px_allshots(:,:,sh) = POCS_projection(x(:,:,sh), sense_map, kspa(:,:,:,sh), trj);
    end
    
    %shot-average
    fpx_allshots = fft2d(px_allshots);
    fpx_allshots_lf = bsxfun(@times, fpx_allshots, W);
    px_allshots_lr = ifft2d(fpx_allshots_lf);
    
    %remove lowresolution phase, i.e. make shots have the same lowres. phase
%     px_aver = mean(px_allshots,3);
    px_aver = mean(px_allshots.* (conj(px_allshots_lr) ./ abs(px_allshots_lr)),3); 
    
    %update extimation
    x_esti_old = x_esti;
    x_esti = x_esti_old + pars.POCS.lamda.*(px_aver - x_esti_old);      
    
    %update initial point for next step (recover true shot by shot phase)
    x = bsxfun(@times, x_esti, (px_allshots_lr./abs(px_allshots_lr)));
%     x = repmat(x_esti,[1 1 42]);
    
    if(it>1)
        residual_error = trace((x_esti-x_esti_old)'*(x_esti-x_esti_old)) / trace(x_esti_old'*x_esti_old);
    end

%     disp(['nit: ', num2str(it)]); 
%     disp(['residual_error: ', num2str(residual_error)]); 
    
    figure(11); 
    subplot(221); imshow(abs(x_esti),[]); title('esti. mag.');
    subplot(222); imshow(angle(x_esti),[-pi pi]); title('esti. phase');drawnow(); 
    subplot(212); plot(it, log10(residual_error),'ro'); xlim([1 pars.POCS.nit]); ylim([-11 0]); title('log10-residual error'); hold on;
    
    it = it +1;
    
    
end
cla

sol =  x_esti;

end