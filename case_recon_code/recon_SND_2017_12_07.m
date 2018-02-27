
clear; clc; close all
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_12_07_3D_spiral')
%% trajectory calculation (Opt 1)-- OFFCENTER MEASURE
close all; clear; clc;
trj_save_fn = 'traj_for_cor.mat';
trajectory_measure_distance = 15; %in mm
spira_3D_trjectory_calculation_normal_6steps(trj_save_fn, trajectory_measure_distance);
disp('-finished- ');

%% trajectory calculation (Opt 2)-- IDEAL WAVEFORM INTERPOLATION
trj_save_fn = 'traj_no_3.mat';

gr_waveform = dlmread('no3_external_spiral_GRwaveforms_FOV200_120_vox_8_6.dat');  %in mT/m

k_all = gr_waveform .* 0;
gr_dwell    = 6.4e-6;   %in s
gamma       = 4257.6;   %Hz/G
G_cm2mT_m   = 10; %unit convert from G/cm to mT/m
for idx = 2:length(k_all)
    for dir = 1:3
        k_all(idx,dir) =  k_all(idx-1,dir) + gamma .* gr_dwell .* (gr_waveform(idx-1,dir)./ G_cm2mT_m);
    end
end
trj_meas_kx_ori = k_all(:,1);
trj_meas_ky_ori = k_all(:,2);
trj_meas_kz_ori = k_all(:,3);

figure(22);
subplot(211); plot(gr_waveform); title('gradient waveform'); ylabel('mT/m')
subplot(212); plot(k_all); title('k space trajectory'); ylabel('1/cm')
save(trj_save_fn, 'trj_meas_kx_ori', 'trj_meas_ky_ori', 'trj_meas_kz_ori' );
%% SET path for all the following steps
clear; close all; clc

data_fn = 'ca_07122017_1822173_9_2_wip_3d_cardiacV4.raw';
sense_ref_fn = 'ca_07122017_1807158_1000_13_wip_senserefscanV4.raw';
coil_survey_fn  = 'ca_07122017_1756277_1000_2_wip_coilsurveyscanV4.raw';

trj_mat_fn = 'traj_no_3.mat';

%% Spiral Nav. data loading
disp('spiral data loading...')
virtual_coil = 0; % 0 for no coil compression
[nav_k_spa_data, Nav_VirtualCoilMartix] = spiral_kspa_data_read(data_fn, virtual_coil);

disp('-finished- ');
%% Interpolate trajecotry
AQ_samples = 3181;
AQ_interval = 0.00702;

load(trj_mat_fn);
trj_meas_kx = interp1q([0:length(trj_meas_kx_ori)-1]'.*0.0064,trj_meas_kx_ori,[0:AQ_samples-1]'.*AQ_interval);
trj_meas_ky = interp1q([0:length(trj_meas_ky_ori)-1]'.*0.0064,trj_meas_ky_ori,[0:AQ_samples-1]'.*AQ_interval);
trj_meas_kz = interp1q([0:length(trj_meas_kz_ori)-1]'.*0.0064,trj_meas_kz_ori,[0:AQ_samples-1]'.*AQ_interval);
assert(sum(isnan(trj_meas_kx))==0); assert(sum(isnan(trj_meas_ky))==0); assert(sum(isnan(trj_meas_kz))==0);

disp(['gr dur: ', num2str((length(trj_meas_kx_ori)-1).*0.0064),'; AQ dur: ', num2str((AQ_samples-1).*AQ_interval)]);

figure(23);
subplot(311);
plot([0:length(trj_meas_kx_ori)-1]'.*0.0064,trj_meas_kx_ori,'b--' ); 
hold on; plot([0:AQ_samples-1]'.*AQ_interval,trj_meas_kx,'r--' ); 
subplot(312);
plot([0:length(trj_meas_ky_ori)-1]'.*0.0064,trj_meas_ky_ori,'b--' ); 
hold on; plot([0:AQ_samples-1]'.*AQ_interval,trj_meas_ky,'r--' ); 
subplot(313);
plot([0:length(trj_meas_kz_ori)-1]'.*0.0064,trj_meas_kz_ori,'b--' ); 
hold on; plot([0:AQ_samples-1]'.*AQ_interval,trj_meas_kz,'r--' ); 

save(trj_mat_fn, 'trj_meas_kz', 'trj_meas_ky', 'trj_meas_kx' ,'-append');
%% Spiral NUFFT recon.
disp(' Spiral NUFFT recon...');
save_mat_fn = 'Sc9.mat';
close all;
[kx_length ch_nr shot_nr, dyn_nr] = size(nav_k_spa_data);

offcenter_xy = [0 0]; 
FOV_xy = [200 200];
% nav_im_recon_nufft = [];
dyn_recon = 1:dyn_nr;
for d = 1:length(dyn_recon)
    tic
    dyn  = dyn_recon(d);
    disp(['dynamic: ',num2str(dyn)]);
    %=============== recon parameters =========================
    recon_par.ignore_kz = 0;
    recon_par.acq_dim = [25 25 20];  
    recon_par.recon_dim  = [25 25 20];
    recon_par.dyn_nr = dyn;
    recon_par.skip_point = 100 ;
    recon_par.end_point = []; %or []: till the end;
    recon_par.interations = 5;
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
        nav_sense_Psi = [];
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
    nav_im_recon_nufft_1dyn = NUFFT_3D_recon(nav_k_spa_data(:,:,:,:),trj_mat_fn,recon_par, nav_sense_map(:,:,:,:), nav_sense_Psi,offcenter_xy, FOV_xy);
    nav_im_recon_nufft(:,:,:,:,:,dyn) = nav_im_recon_nufft_1dyn;
    save(save_mat_fn, 'nav_im_recon_nufft','-append'); 
    
    
    elaps_t=toc;
    msg = sprintf(['SoSNav recon finishted for {', data_fn,'} ; ...dynamic %d ; duration %f; s', 10, 'Saved as ', save_mat_fn],d, elaps_t);
%     sendmail_from_yahoo('q.zhang@amc.nl','Matlab Message',msg);
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
