
clear; clc; close all
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_11_27_3D_spiral')
%% trajectory calculation
close all; clear; clc;
trj_save_fn = 'traj_for_Sc6.mat';
trajectory_measure_distance = 15; %in mm
spira_3D_trjectory_calculation_normal_6steps(trj_save_fn, trajectory_measure_distance);
disp('-finished- ');

%% SET path for all the following steps
clear; close all; clc

data_fn = 'sp_27112017_1643025_6_2_wip_3d_brain_cardiacV4.raw';
sense_ref_fn = 'sp_27112017_1642355_1000_28_wip_senserefscanV4.raw';
coil_survey_fn  = 'sp_27112017_1642147_1000_25_wip_coilsurveyscanV4.raw';

trj_mat_fn = 'traj_for_Sc6.mat';

%% Spiral Nav. data loading
disp('spiral data loading...')
virtual_coil = 0; % 0 for no coil compression
[nav_k_spa_data, Nav_VirtualCoilMartix] = spiral_kspa_data_read(data_fn, virtual_coil);

disp('-finished- ');
%% Spiral NUFFT recon.
disp(' Spiral NUFFT recon...');
save_mat_fn = 'Sc6.mat';
close all;
[kx_length ch_nr shot_nr, dyn_nr] = size(nav_k_spa_data);

offcenter_xy = [0 0]; 
FOV_xy = [300 300];
% nav_im_recon_nufft = [];
dyn_recon = 1:dyn_nr;
for d = 1:length(dyn_recon)
    tic
    dyn  = dyn_recon(d);
    disp(['dynamic: ',num2str(dyn)]);
    %=============== recon parameters =========================
    recon_par.ignore_kz = 0;
    recon_par.acq_dim = [30 30 15];  
    recon_par.recon_dim  = [30 30 15];
    recon_par.dyn_nr = dyn;
    recon_par.skip_point = 108 ;
    recon_par.end_point = []; %or []: till the end;
    recon_par.interations = 10;
    recon_par.lamda = 0.1;
    recon_par.recon_all_shot = 1;
    recon_par.sense_map_recon =1; 
    recon_par.update_SENSE_map = 0;
    recon_par.sense_calc_method = 'external'; %'ecalib' or 'external'
    recon_par.sense_os = [1 FOV_xy(1)/FOV_xy(2)];  %oversampling in x and y: control sense FOV
    recon_par.data_fn = data_fn;
    recon_par.sense_ref = sense_ref_fn;
    recon_par.coil_survey = coil_survey_fn;
    
    recon_par.channel_by_channel = 1;
    recon_par.channel_by_channel = recon_par.channel_by_channel .* (1-recon_par.sense_map_recon );
    %========================  END  =========================
     if(~exist('nav_sense_map', 'var')&&recon_par.sense_map_recon)
        recon_par.update_SENSE_map = 1;
     end
    
    if(recon_par.update_SENSE_map)
        [nav_sense_map, nav_sense_Psi] = calc_sense_map(recon_par.data_fn, recon_par.sense_ref,  recon_par.coil_survey, recon_par.recon_dim,recon_par.sense_calc_method, recon_par.sense_os);
        %compress sense map and sense_Psi
        if(exist('Nav_VirtualCoilMartix','var'))
            if(~isempty(Nav_VirtualCoilMartix))
                [nav_sense_map, nav_sense_Psi] = compress_sense_map_Psi(Nav_VirtualCoilMartix, nav_sense_map, nav_sense_Psi);
            end
        end
    end
    
    if(recon_par.sense_map_recon == 0)
        nav_sense_map = ones([recon_par.recon_dim ch_nr]);
    end
    nav_sense_map = normalize_sense_map(nav_sense_map);
    
    
    
    if(~exist('nav_sense_Psi','var'))
        nav_sense_Psi = [];
    end
    nav_im_recon_nufft_1dyn = NUFFT_3D_recon(nav_k_spa_data,trj_mat_fn,recon_par, nav_sense_map, nav_sense_Psi,offcenter_xy, FOV_xy);
    nav_im_recon_nufft(:,:,:,:,:,dyn) = nav_im_recon_nufft_1dyn;
    save(save_mat_fn, 'nav_im_recon_nufft','-append'); 
    
    
    elaps_t=toc;
    msg = sprintf(['SoSNav recon finishted for {', data_fn,'} ; ...dynamic %d ; duration %f; s', 10, 'Saved as ', save_mat_fn],d, elaps_t);
    sendmail_from_yahoo('q.zhang@amc.nl','Matlab Message',msg);
end
% nav_sense_map = circshift(nav_sense_map, round(17.26/115.00*size(nav_sense_map,1)));
% nav_im_recon_nufft = circshift(nav_im_recon_nufft, -1*round(17.26/115.00*size(nav_sense_map,1)));
dyn = 1;
figure(801); immontage4D(angle(squeeze(nav_im_recon_nufft(:,:,:,:,:,dyn))),[-pi pi]); colormap jet; 
figure(802); immontage4D(abs(squeeze(nav_im_recon_nufft(:,:,:,:,:,dyn))),[]); 
phase_diff = angle(squeeze(bsxfun(@times,  nav_im_recon_nufft, exp(-1i*angle(nav_im_recon_nufft(:,:,:,1,1,:))))));
figure(803); immontage4D(squeeze(phase_diff(:,:,:,:,dyn)),[-pi pi]); colormap jet;

if(recon_par.channel_by_channel)
    nav_im_ch_by_ch = nav_im_recon_nufft_1dyn;
end

if(exist('nav_sense_map','var')&&exist('nav_im_ch_by_ch','var'))
    figure(804); 
    slice1=ceil(size(nav_im_ch_by_ch,3)/2);
    slice2=ceil(size(nav_sense_map,3)/2);
    subplot(121); montage(abs(nav_im_ch_by_ch(:,:,slice1,:)),'displayrange',[]); title('Check if they are match!'); xlabel('channel-by-channel');
    subplot(122); montage(abs(nav_sense_map(:,:,slice2,:)),'displayrange',[]); xlabel('sense');
end


disp('-finished- ');
