function nav_im_recon_nufft = NUFFT_3D_recon(data_fn,trj_fn,recon_par)


load(trj_fn); load(data_fn);

if(isempty(recon_par.end_point))
    recon_par.end_point = length(trj_meas_kx);
end

selected_point = recon_par.skip_point+1:recon_par.end_point;

%calc sens_maps
if(recon_par.sense_map_recon)
    %get k space
    if(exist('sens_map')&(~recon_par.update_SENSE_map))
    else
        if(strcmp(recon_par.sense_calc_method, 'ecalib'))
            
            sens_map = get_sense_map_ecalib(recon_par.data_fn,recon_par.recon_dim );
            
        elseif (strcmp(recon_par.sense_calc_method, 'external'))
            
            sens_map = get_sense_map_external(recon_par.sense_ref, recon_par.data_fn, recon_par.coil_survey, recon_par.recon_dim);
        end
        save(data_fn, 'sens_map','-append');
    end
    
else
    sens_map = ones(recon_par.recon_dim);
end




trj_meas_kx_t = squeeze(trj_meas_kx(selected_point,1));
trj_meas_ky_t = squeeze(trj_meas_ky(selected_point,1));
trj_meas_kz_t = squeeze(trj_meas_kz(selected_point,1));

clear trj_meas_kx trj_meas_ky trj_meas_kz k_spa_data
trj_meas_kx = trj_meas_kx_t';
trj_meas_ky = trj_meas_ky_t';
trj_meas_kz = trj_meas_kz_t';

clear trj_meas_kx_t trj_meas_ky_t trj_meas_kz_t sig_bart_t

scale_foctor_xy = max(2*pi/(max(trj_meas_kx)-min(trj_meas_kx)),2*pi/(max(trj_meas_ky)-min(trj_meas_ky)));
trj_meas_kx_scaled = trj_meas_kx * scale_foctor_xy;
trj_meas_ky_scaled = trj_meas_ky * scale_foctor_xy;

scale_foctor_z = 2*pi/(max(trj_meas_kz)-min(trj_meas_kz));

if(recon_par.ignore_kz == 1)
    scale_foctor_z = 0;
end
trj_meas_kz_scaled = trj_meas_kz * scale_foctor_z;


figure(8);
plot(trj_meas_kx_scaled); hold on; plot(trj_meas_ky_scaled,'r');  plot(trj_meas_kz_scaled,'k');legend('measured kx (a.u.)','measured ky (a.u.)','measured kz (a.u.)');

figure(801);
plot3(trj_meas_kx_scaled, trj_meas_ky_scaled, trj_meas_kz_scaled); title('-pi to pi')


trj_nufft = double(cat(1, trj_meas_kx_scaled, trj_meas_ky_scaled, trj_meas_kz_scaled));
trj_nufft = trj_nufft';
clear im_recon_nufft nav_im_recon_nufft



% for shot_nr = 1: size(nav_k_spa_data,3)
for shot_nr = 1:1
    
    shot_nr
    
    sig_kspa = nav_k_spa_data(selected_point,:,shot_nr,recon_par.dyn_nr);
    
    if(recon_par.sense_map_recon) %all in one recon
        sig_nufft = col(double(sig_kspa));
        A=nuFTOperator(trj_nufft,recon_par.recon_dim,sens_map,6);
        
        % simple inverse
        nav_im_recon_nufft = A'*sig_nufft;
        figure(702); montage(permute(abs(nav_im_recon_nufft), [1 2 4 3]),'displayrange',[]); title('direct inverse');
        
        
        
        %call CG-SENSE with L2-norm regularization
             nav_im_recon_nufft(:,:,:,:,shot_nr)=regularizedReconstruction(A,sig_nufft,@L2Norm,0.5,'maxit',recon_par.interations);
%         nav_im_recon_nufft(:,:,:,:,shot_nr)=regularizedReconstruction(A,sig_nufft,'maxit',recon_par.interations);
        
    else %ch by ch recon
        
        for ch =1:size(nav_k_spa_data,2)
            sig_nufft = double(sig_kspa(:, ch));
            sig_nufft = sig_nufft';
            A=nuFTOperator(trj_nufft,recon_par.recon_dim,sens_map,6);
            
            % simple inverse
            %                 nav_im_recon_nufft(:,:,:,ch) = A'*sig_nufft';
            
            
            
            %call CG-SENSE with L2-norm regularization
                 nav_im_recon_nufft(:,:,:,ch)=regularizedReconstruction(A,sig_nufft',@L2Norm,0.5,'maxit',recon_par.interations);
%             nav_im_recon_nufft(:,:,:,ch,shot_nr)=regularizedReconstruction(A,sig_nufft','maxit',recon_par.interations);
        end
        % im_recon_nufft = flipdim(flipdim(im_recon_nufft,1),2);
    end
end


display_shot_nr = 1;

if(recon_par.sense_map_recon) %all in one recon
    figure(37); montage(permute(abs(nav_im_recon_nufft), [1 2 4 3]),'displayrange',[]); title('all slices rss'); xlabel('P'); ylabel('M')
    figure(38); montage(permute(abs(nav_im_recon_nufft), [2 3 4 1]),'displayrange',[]); title('all slices rss'); xlabel('S'); ylabel('P')
else
    
    figure(35); montage(permute(abs(squeeze(nav_im_recon_nufft(:,:,:,6,display_shot_nr))), [1 2 4 3]),'displayrange',[]); title('all slices ch6')
    figure(36); montage(-1* permute(angle(squeeze(nav_im_recon_nufft(:,:,:,6,display_shot_nr))), [1 2 4 3]),'displayrange',[-pi pi]); colormap jet; title('all slices ch6')
    
    nav_im_recon_nufft_rss = squeeze(bart('rss 8', nav_im_recon_nufft));
    figure(37); montage(permute(abs(nav_im_recon_nufft_rss), [1 2 4 3]),'displayrange',[]); title('all slices rss'); xlabel('P'); ylabel('M')
    figure(38); montage(permute(abs(nav_im_recon_nufft_rss), [2 3 4 1]),'displayrange',[]); title('all slices rss'); xlabel('S'); ylabel('P')
end
end