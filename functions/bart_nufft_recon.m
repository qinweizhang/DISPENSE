function [reco2, igrid, igrid_rss] = bart_nufft_recon(data_fn, trj_fn, bart_recon)


load(data_fn); load(trj_fn);


if(isempty(bart_recon.end_point))
    bart_recon.end_point = length(trj_meas_kx);
end


selected_point = [bart_recon.skip_point+1:bart_recon.end_point];
trj_meas_kx_t = trj_meas_kx(selected_point,bart_recon.diffusion_nr);
trj_meas_ky_t = trj_meas_ky(selected_point,bart_recon.diffusion_nr);
trj_meas_kz_t = trj_meas_kz(selected_point,bart_recon.diffusion_nr);

clear trj_meas_kx trj_meas_ky trj_meas_kz
trj_meas_kx = trj_meas_kx_t;
trj_meas_ky = trj_meas_ky_t;
trj_meas_kz = trj_meas_kz_t;


scale_foctor_kx_ky = max([bart_recon.trj_scale_dim(1)/(max(trj_meas_kx)-min(trj_meas_kx)),...
    bart_recon.trj_scale_dim(2)/(max(trj_meas_ky)-min(trj_meas_ky))]);
scale_foctor_kz = bart_recon.trj_scale_dim(3)/(max(trj_meas_kz)-min(trj_meas_kz));

trj_meas_kx_scaled = trj_meas_kx * scale_foctor_kx_ky;
trj_meas_ky_scaled = trj_meas_ky * scale_foctor_kx_ky;
trj_meas_kz_scaled = trj_meas_kz * scale_foctor_kz;


figure(8);
subplot(121); plot(trj_meas_kx_scaled); hold on; plot(trj_meas_ky_scaled,'r');plot(trj_meas_kz_scaled,'k'); legend('measured kx (a.u.)','measured ky (a.u.)','measured kz (a.u.)');
subplot(122); plot3(trj_meas_kx_scaled, trj_meas_ky_scaled, trj_meas_kz_scaled); xlabel('measured kx (a.u.)'); ylabel('measured ky (a.u.)'); zlabel('measured kz (a.u.)')

clear igrid sig_bart

if(bart_recon.ignor_kz)
    trj_bart = double(cat(2, trj_meas_kx_scaled, trj_meas_ky_scaled, zeros(size(trj_meas_kx_scaled))));
else
    trj_bart = double(cat(2, trj_meas_kx_scaled, trj_meas_ky_scaled, trj_meas_kz_scaled));
end

sig_bart(1,:,1,:) = nav_k_spa_data(selected_point,:,bart_recon.nsa_nr,bart_recon.shot_nr,bart_recon.diffusion_nr);

bart_command_nufft = sprintf( 'nufft -i -d%d:%d:%d -l0.01',bart_recon.recon_dim(1), bart_recon.recon_dim(2), bart_recon.recon_dim(3))
igrid=((bart(bart_command_nufft,trj_bart',sig_bart))); % nufft take trj  in [3,kx_range,[ky_range],[kz_range]]; kspa_data in [1,kx_range,[ky_range],[kz_range],ch]

igrid_rss = bart('rss 8',igrid);
slice_id = 3;
figure(21); montage(abs(igrid(:,:,slice_id,:)),'Displayrange',[],'size',[4 4]); title('mag. igrid: slice 5, all channels')
figure(22); montage(angle(igrid(:,:,slice_id,:)),'Displayrange',[-pi pi],'size',[4 4]); colormap 'jet'; title('phase. igrid: : slice 5, all channel')
figure(23); montage(abs(permute(igrid_rss,[1 2 4 3])),'Displayrange',[]); title('mag. igrid_rss')

%% %-----BART PICS-------%
if(bart_recon.PICS)
    
    % get sens maps
    if(exist('nav_bart_sens_map')&(~bart_recon.update_SENSE_map))
    else
        if(strcmp(bart_recon.sense_calc_method, 'ecalib'))
            
            % reconstruct low-resolution image and transform back to k-space
            bart_command_1 = sprintf('nufft -i -d%d:%d:%d -l0.01',24, 24, 10); %now fixed
            lowres_img = bart(bart_command_1, trj_bart',sig_bart );
            lowres_ksp = bart('fft -u 7', lowres_img);
            
            % figure(202); montage(abs(lowres_img(:,:,2,:)),'displayrange',[]);
            
            % zeropad to full size
            clear fullres_ksp
            bart_command_2 = sprintf('resize -c 0 %d 1 %d 2 %d',bart_recon.recon_dim(1),bart_recon.recon_dim(2),bart_recon.recon_dim(3))
            ksp_zerop = bart(bart_command_2, lowres_ksp);
            fullres_ksp(:,:,:,:) = bart('fft -u 7', igrid);
            
            % ESPIRiT calibration
            clear nav_bart_sens_map
            % sens = bart('ecalib -m1', ksp_zerop); %low res sense
            nav_bart_sens_map = bart('ecalib -m1', fullres_ksp); %high res sense
            
        elseif (strcmp(bart_recon.sense_calc_method, 'external'))
            
            nav_bart_sens_map = get_sense_map_external(bart_recon.sense_ref, bart_recon.data_fn, bart_recon.coil_survey, bart_recon.recon_dim);
        end
        
        
        save(data_fn, 'nav_bart_sens_map','-append');
    end
    
    %PICS
    bart_command_3 = sprintf('pics -S -r0.001 -t')
    reco2 = bart(bart_command_3, trj_bart', sig_bart, nav_bart_sens_map);
    
    
end


end