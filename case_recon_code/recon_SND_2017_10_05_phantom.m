
clear; clc; close all
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_10_05_SND_phantom')
%% trajectory calculation
close all; clear; clc;
trj_save_fn = 'traj_for_Sc19.mat';
trajectory_measure_distance = 15; %in mm
spira_3D_trjectory_calculation(trj_save_fn, trajectory_measure_distance);
disp('-finished- ');

%% SET path for all the following steps
clear; close all; clc

data_fn = 'sn_05102017_1741096_5_2_wip_sc17_dpsti_sosad_linearV4.raw';
sense_ref_fn = 'sn_05102017_1720266_1000_7_wip_senserefscanV4.raw';
coil_survey_fn  = 'sn_05102017_1717314_1000_2_wip_coilsurveyscanV4.raw';

trj_mat_fn = 'traj_Sc25_26_27_for_Sc17.mat';



%% Spiral Nav. data loading
disp('spiral Nav. data loading...')
nav_k_spa_data = nav_kspa_data_read(data_fn);

disp('-finished- ');
%% Spiral NUFFT recon.
disp(' Spiral NUFFT recon...');
close all;
[kx_length ch_nr shot_nr, dyn_nr] = size(nav_k_spa_data);

nav_im_recon_nufft = [];
for dyn = 1:dyn_nr
    dyn
    %=============== recon parameters =========================
    recon_par.ignore_kz = 0;
    recon_par.acq_dim = [26 26 10];  
    recon_par.recon_dim  = [26 26 10];
    recon_par.dyn_nr = dyn;
    recon_par.skip_point = 0 ;
    recon_par.end_point = []; %or []: till the end;
    recon_par.interations = 5;
    recon_par.lamda = 0.5;
    recon_par.recon_all_shot = 1;
    recon_par.sense_map_recon =1; 
    recon_par.update_SENSE_map = 0;
    recon_par.sense_calc_method = 'external'; %'ecalib' or 'external'
    recon_par.sense_os = [1 1];  %oversampling in x and y: control sense FOV
    recon_par.data_fn = data_fn;
    recon_par.sense_ref = sense_ref_fn;
    recon_par.coil_survey = coil_survey_fn;
    
    recon_par.channel_by_channel = 1;
    recon_par.channel_by_channel = recon_par.channel_by_channel .* (1-recon_par.sense_map_recon );
    %========================  END  =========================
     if(~exist('nav_sense_map', 'var')&&recon_par.sense_map_recon)
        recon_par.update_SENSE_map = 1;
     end
    
    nav_sense_Psi = [];
    if(recon_par.update_SENSE_map)
        [nav_sense_map, nav_sense_Psi] = calc_sense_map(recon_par.data_fn, recon_par.sense_ref,  recon_par.coil_survey, recon_par.recon_dim,recon_par.sense_calc_method, recon_par.sense_os);
    end
    
    if(recon_par.sense_map_recon == 0)
        nav_sense_map = ones([recon_par.recon_dim ch_nr]);
    end
    nav_sense_map = normalize_sense_map(nav_sense_map);
    
    offcenter_xy = [0 0]; 
    FOV_xy = [210 210];
    
    nav_im_recon_nufft_1dyn = NUFFT_3D_recon(nav_k_spa_data,trj_mat_fn,recon_par, (nav_sense_map));
    nav_im_recon_nufft = cat(6, nav_im_recon_nufft, nav_im_recon_nufft_1dyn);
end
% nav_sense_map = circshift(nav_sense_map, round(17.26/115.00*size(nav_sense_map,1)));
% nav_im_recon_nufft = circshift(nav_im_recon_nufft, -1*round(17.26/115.00*size(nav_sense_map,1)));
figure(801); immontage4D(angle(squeeze(nav_im_recon_nufft)),[-pi pi]); colormap jet; 
figure(8011); immontage4D(angle(permute(squeeze(nav_im_recon_nufft),[1 3 2 4])), [-pi pi]); colormap jet; 

figure(802); immontage4D(abs(squeeze(nav_im_recon_nufft)),[]); 

phase_diff = angle(squeeze(bsxfun(@times,  nav_im_recon_nufft, exp(-1i*angle(nav_im_recon_nufft(:,:,:,:,1))))));
figure(803); immontage4D(squeeze(phase_diff),[-pi pi]); colormap jet;
figure(804); immontage4D(permute(squeeze(phase_diff),[1 3 2 4]),[-pi pi]); colormap jet;


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

%% 3D sprial navigator data phase unwrapping
current_mat_file = 'Sc5.mat'
nav_im_recon_nufft=permute(nav_im_recon_nufft,[1 3 2 4]);

nav_im_recon_nufft_diff = bsxfun(@rdivide, nav_im_recon_nufft,  exp(1i.*angle(nav_im_recon_nufft(:,:,:,1))));
nav_im_recon_nufft_diff(find(isinf(nav_im_recon_nufft_diff)))=0;
nav_im_recon_nufft_diff(find(isnan(nav_im_recon_nufft_diff)))=0;


[nav_dimx, nav_dimy, nav_dimz, nav_dyn] = size(nav_im_recon_nufft_diff);



%------------------------------------------------------PHSE UNWRAPPING----------------------------------%

Manually_seed_point = 0;
for dyn = 1:nav_dyn
    for z = 1:nav_dimz
        nav_ima_phase_unwrapped(:,:,z,dyn) = GoldsteinUnwrap2D_r1_func(squeeze(nav_im_recon_nufft_diff(:,:,z,dyn)), Manually_seed_point);
        figure(22);
        subplot(121); imshow(squeeze(angle(nav_im_recon_nufft_diff(:,:,z,dyn))),[-pi pi]); title('L: wraped R: unwrapped ') ;  colorbar
        subplot(122); imshow(squeeze(nav_ima_phase_unwrapped(:,:,z,dyn)),[-pi pi]); colormap jet;title(['dyn: ',num2str(dyn),'z: ',num2str(z)]);  colorbar
        drawnow;
        
    end
end


% ===================display wrapped phase
phase_wrapped = squeeze(angle(nav_im_recon_nufft_diff));


figure(601);
subplot(121); montage(permute(squeeze(phase_wrapped(:,:,round(nav_dimz/2),:)),[1 2 4 3]),'displayrange',[-pi pi]); colormap jet;
subplot(122); montage(permute(squeeze(nav_ima_phase_unwrapped(:,:,round(nav_dimz/2),:)),[1 2 4 3]),'displayrange',[-pi pi]); colormap jet;
% ===================display wrapped phase difference
phase_test_ch = phase_wrapped(:,:,round(nav_dimz/2),:);
for shot_idx = 1:27
    phase_test_ch_diff(:,:,1,shot_idx) = phase_test_ch(:,:,1,shot_idx) -phase_test_ch(:,:,1,4) ;
end
figure(602); montage(phase_test_ch_diff,'displayrange',[-pi pi]); colormap jet; title(['all shots phase differece from ch', num2str(dyn)]);


% ====================display unwrapped phase difference
phase_test_ch = nav_ima_phase_unwrapped(:,:,round(nav_dimz/2),:);
for shot_idx = 1:27
    phase_test_ch_diff(:,:,1,shot_idx) = phase_test_ch(:,:,1,shot_idx) -phase_test_ch(:,:,1,4) ;
end
figure(604); montage(phase_test_ch_diff,'displayrange',[-5*pi 5*pi]); colormap jet; title(['unwrapped phase difference:all shots phase differece from ch', num2str(dyn)]);



nav_im_recon_nufft = permute(squeeze(nav_im_recon_nufft),[1 2 3 4]);
nav_ima_phase_unwrapped =  permute(squeeze(nav_ima_phase_unwrapped),[1 2 3 4]);
% nav_im_recon_nufft =  nav_im_recon_nufft(:,:,6:7,:);
% nav_ima_phase_unwrapped =  nav_ima_phase_unwrapped(:,:,6:7,:);
save(current_mat_file,'nav_ima_phase_unwrapped', 'nav_im_recon_nufft');

%% --------------Navigator processing--------------------
clc; close all;

nav_kspa_recon_nufft = bart('fft 3',nav_im_recon_nufft); 
[kx, ky, ch, shot, dyn] = size(nav_kspa_recon_nufft);
save(current_mat_file,'nav_kspa_recon_nufft','-append');

[ nav_ima_phase_unwrapped_diff,fitted_nav_ima_phase, linear_phase_xy,global_phase] =  easy_rigidMotion_parameter_calculation( nav_ima_phase_unwrapped, nav_im_recon_nufft);
save(current_mat_file, 'nav_ima_phase_unwrapped_diff','fitted_nav_ima_phase', 'linear_phase_xy','global_phase','-append');

%%
edit analyze_input_phase_error
















%%  

% ==========================================================================================
% ===                                                                                 ======
% ===                                                                                 ======
% ===                                                                                 ======
% ===                                                                                 ======
% ===             TSE        BELOW                                                    ======
% ===                                                                                 ======
% ===                                                                                 ======
% ===                                                                                 ======
% ===                                                                                 ======
% ===                                                                                 ======
% ==========================================================================================


%% TSE data sorting and default recon
close all; clc
disp(' TSE data sorting and default recon...')

parameter2read.dyn = [];

[ima_k_spa_data,TSE.ky_matched,TSE.kz_matched,TSE.shot_matched, TSE.ch_dim,ima_kspa_sorted, ima_default_recon, TSE_sense_map, TSE.kxrange, TSE.kyrange, TSE.kzrange] = ...
    TSE_data_sortting(data_fn, sense_ref_fn, coil_survey_fn,parameter2read);

TSE

figure(606); immontage4D(permute(abs(ima_default_recon(:,:,:,:)),[1 2 4 3]), [100 1500]);
disp('-finished- ');

%% TSE data non-rigid phase error correction (iterative) CG_SENSE

TSE.kxrange = [-224 -1]; %consider now the ima_k_spa_data is oversampled in kx
TSE.kyrange = [-112 -1]; %consider now the ima_k_spa_data is oversampled in kx

TSE.dyn_dim = 27;

%parameters for DPsti_TSE_phase_error_cor
pars.sense_map = 'external';  % external or ecalib 

pars.data_fn = data_fn;
pars.sense_ref = sense_ref_fn;
pars.coil_survey = coil_survey_fn;
pars.enabled_ch = [4]; % take the 4th channel
for dyn =1:27
pars.b0_shots = [1:18] + (dyn-1)*18;


TSE.kx_dim = TSE.kxrange(2) - TSE.kxrange(1) + 1;
TSE.ky_dim = TSE.kyrange(2) - TSE.kyrange(1) + 1;  %max_ky * 2 + 1;
TSE.kz_dim = TSE.kzrange(2) - TSE.kzrange(1) + 1;

% Preprocessing on kspace data for b=0
shots_per_dyn =  max(TSE.shot_matched)/TSE.dyn_dim;
disp('Preprocessing on kspace data for b=0');
if(isempty(pars.b0_shots))
    b0_shots_range = 1:shots_per_dyn; %by default, the first dynamic
else
    b0_shots_range = pars.b0_shots;
end

tic;

kspa_b0 = sort_k_spa_sh_by_sh(ima_k_spa_data, b0_shots_range, TSE, pars);

kspa_b0_combshot = sum(kspa_b0, 5)./sum(abs(kspa_b0)>0, 5); %4D b0 kspace [kx ky kz nc]; non-zeros average
kspa_b0_combshot(find(isnan(kspa_b0_combshot)))=0; kspa_b0_combshot(find(isinf(kspa_b0_combshot)))=0;

im_b0_ch_by_ch=ifft3d(kspa_b0_combshot);
figure(1); immontage4D(abs(im_b0_ch_by_ch),[]);
figure(2); immontage4D(angle(im_b0_ch_by_ch),[-pi pi]);colormap jet

TSE_all_temp(:,:,:,dyn) = im_b0_ch_by_ch(:,:,:); 
end

TSE_all=permute(squeeze(TSE_all_temp(:,:,14:15,:)), [1 2 3 4]);
figure(3); immontage4D(angle(TSE_all),[-pi pi]);colormap jet


toc;

%% 3D sprial TSE data phase unwrapping

current_mat_file = 'Sc23.mat'

TSE_im_recon_nufft_diff = bsxfun(@rdivide, TSE_all, exp(1i.*angle(TSE_all(:,:,:,1))));
TSE_im_recon_nufft_diff(find(isinf(TSE_im_recon_nufft_diff)))=0;
TSE_im_recon_nufft_diff(find(isnan(TSE_im_recon_nufft_diff)))=0;
figure(4); immontage4D(angle(TSE_im_recon_nufft_diff),[-pi pi]);colormap jet

[nav_dimx, nav_dimy, nav_dimz, nav_dyn] = size(TSE_im_recon_nufft_diff);



%------------------------------------------------------PHSE UNWRAPPING----------------------------------%

Manually_seed_point = 0;
for dyn = 1:nav_dyn
    for z = 1:nav_dimz
        TSE_ima_phase_unwrapped(:,:,z,dyn) = GoldsteinUnwrap2D_r1_func(squeeze(TSE_im_recon_nufft_diff(:,:,z,dyn)), Manually_seed_point);
        figure(22);
        subplot(121); imshow(squeeze(angle(TSE_im_recon_nufft_diff(:,:,z,dyn))),[-pi pi]); title('L: wraped R: unwrapped ') ;  colorbar
        subplot(122); imshow(squeeze(TSE_ima_phase_unwrapped(:,:,z,dyn)),[-pi pi]); colormap jet;title(['dyn: ',num2str(dyn),'z: ',num2str(z)]);  colorbar
        drawnow;
        
    end
end


% ===================display wrapped phase
phase_wrapped = angle(TSE_all);

figure(601);
subplot(121); montage(permute(squeeze(phase_wrapped(:,:,round(nav_dimz/2),:)),[1 2 4 3]),'displayrange',[-pi pi]); colormap jet;
subplot(122); montage(permute(squeeze(TSE_ima_phase_unwrapped(:,:,round(nav_dimz/2),:)),[1 2 4 3]),'displayrange',[-pi pi]); colormap jet;
% ===================display wrapped phase difference
phase_test_ch = phase_wrapped(:,:,round(nav_dimz/2),:);
for shot_idx = 1:27
    phase_test_ch_diff(:,:,1,shot_idx) = phase_test_ch(:,:,1,shot_idx) -phase_test_ch(:,:,1,4) ;
end
figure(602); montage(phase_test_ch_diff,'displayrange',[-pi pi]); colormap jet; title(['all shots phase differece from ch', num2str(dyn)]);


% ====================display unwrapped phase difference
phase_test_ch = TSE_ima_phase_unwrapped(:,:,round(nav_dimz/2),:);
for shot_idx = 1:27
    phase_test_ch_diff(:,:,1,shot_idx) = phase_test_ch(:,:,1,shot_idx) -phase_test_ch(:,:,1,4) ;
end
figure(604); montage(phase_test_ch_diff,'displayrange',[-5*pi 5*pi]); colormap jet; title(['unwrapped phase difference:all shots phase differece from ch', num2str(dyn)]);



save(current_mat_file,'TSE_ima_phase_unwrapped', 'TSE_all', '-append'); 

%% --------------TSE phase processing--------------------
clc; close all;

TSE_all_kspa = bart('fft 3',TSE_all); 
[kx, ky, ch,  dyn] = size(TSE_all_kspa);
save(current_mat_file,'TSE_all_kspa','-append');

[ TSE_ima_phase_unwrapped_diff,fitted_TSE_ima_phase, linear_phase_xy_TSE,global_phase_TSE] = easy_rigidMotion_parameter_calculation( TSE_ima_phase_unwrapped, TSE_all);   
save(current_mat_file, 'TSE_ima_phase_unwrapped_diff','fitted_TSE_ima_phase', 'linear_phase_xy_TSE','global_phase_TSE','-append');

%%
linear_phase_xy = linear_phase_xy_TSE;
global_phase = global_phase_TSE;
edit analyze_input_phase_error