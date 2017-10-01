
clear; clc; close all
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_09_24_SND_brain/brain')
%% trajectory calculation
close all; clear; clc;
trj_save_fn = 'traj_for_Sc7.mat';
trajectory_measure_distance = 15; %in mm
spira_3D_trjectory_calculation(trj_save_fn, trajectory_measure_distance);
disp('-finished- ');

%% SET path for all the following steps
clear; close all; clc

data_fn = 'lr_24092017_1551062_4_2_wip_sc3_snd_brain_4b_ppuV4.raw';
sense_ref_fn = 'lr_24092017_1550309_1000_11_wip_senserefscanV4.raw';
coil_survey_fn  = 'lr_24092017_1550042_1000_8_wip_coilsurveyscanV4.raw';

trj_mat_fn = 'traj_for_Sc4.mat';

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
    recon_par.acq_dim = [26 26 11];  
    recon_par.recon_dim  = [26 26 11];
    recon_par.dyn_nr = dyn;
    recon_par.skip_point = 0 ;
    recon_par.end_point = []; %or []: till the end;
    recon_par.interations = 5;
    recon_par.lamda = 2;
    recon_par.recon_all_shot = 1;
    recon_par.sense_map_recon =1; 
    recon_par.update_SENSE_map = 0;
    recon_par.sense_calc_method = 'external'; %'ecalib' or 'external'
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
        nav_sense_map = calc_sense_map(recon_par.data_fn, recon_par.sense_ref,  recon_par.coil_survey, recon_par.recon_dim,recon_par.sense_calc_method);
    end
    
    if(recon_par.sense_map_recon == 0)
        nav_sense_map = ones([recon_par.recon_dim ch_nr]);
    end
    nav_sense_map = normalize_sense_map(nav_sense_map);
    
    
    nav_im_recon_nufft_1dyn = NUFFT_3D_recon(nav_k_spa_data,trj_mat_fn,recon_par, (nav_sense_map));
    nav_im_recon_nufft = cat(6, nav_im_recon_nufft, nav_im_recon_nufft_1dyn);
end
% nav_sense_map = circshift(nav_sense_map, round(17.26/115.00*size(nav_sense_map,1)));
% nav_im_recon_nufft = circshift(nav_im_recon_nufft, -1*round(17.26/115.00*size(nav_sense_map,1)));
figure(801); immontage4D(angle(squeeze(nav_im_recon_nufft(:,:,:,:,:,1))),[-pi pi]); colormap jet; 
figure(802); immontage4D(abs(squeeze(nav_im_recon_nufft(:,:,:,:,:,1))),[]); 
phase_diff = angle(squeeze(bsxfun(@times,  nav_im_recon_nufft, exp(-1i*angle(nav_im_recon_nufft(:,:,:,1,1,:))))));
figure(803); immontage4D(squeeze(phase_diff(:,:,:,:,1)),[-pi pi]); colormap jet;

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
%% -----BART recon -------%
%{ 
disp(' Spiral NUFFT recon (BART)...')
data_fn = 'data_Sc14_3D.mat';
trj_fn = 'traj_Sc25_26_27_for_Sc14.mat';

%======recon parameters ======%
bart_recon.ignor_kz = 0;  % 1for yes. 0 for no
bart_recon.diffusion_nr = 1;
bart_recon.skip_point =0;
bart_recon.end_point = [];
bart_recon.recon_dim = [26 26 10];
bart_recon.trj_scale_dim = [16 16 5];
bart_recon.shot_nr = 1;
bart_recon.nsa_nr = 1;
bart_recon.data_fn = data_fn;
bart_recon.sense_ref = sense_ref_fn;
bart_recon.coil_survey = coil_survey_fn;
bart_recon.update_SENSE_map = 0;
bart_recon.PICS = true;
bart_recon.sense_calc_method = 'external'; %'ecalib' or 'external'
%========end==================%
[reco_pics, igrid, igrid_rss] = bart_nufft_recon(data_fn, trj_fn, bart_recon);
 
figure(204); montage(permute(abs(reco_pics),[1 2 4 3]),'displayrange',[])
figure(205); montage(permute(angle(reco_pics),[1 2 4 3]),'displayrange',[-pi pi]); colormap jet
save(data_fn, 'reco_pics','igrid','igrid_rss','-append');
%}
%% TSE data sorting and default recon
close all; clc
disp(' TSE data sorting and default recon...')

parameter2read.dyn = [];

[ima_k_spa_data,TSE.ky_matched,TSE.kz_matched,TSE.shot_matched, TSE.ch_dim,ima_kspa_sorted, ima_default_recon, TSE_sense_map, TSE.kxrange, TSE.kyrange, TSE.kzrange] = ...
    TSE_data_sortting(data_fn, sense_ref_fn, coil_survey_fn,parameter2read);

TSE

figure(610); immontage4D(permute(abs(ima_default_recon(:,:,:,:)),[1 2 4 3]), [100 500]);
disp('-finished- ');

%% TSE data non-rigid phase error correction (iterative) CG_SENSE
nav_data = reshape(nav_im_recon_nufft, size(nav_im_recon_nufft,1), size(nav_im_recon_nufft, 2), size(nav_im_recon_nufft, 3), max(TSE.shot_matched));

TSE.kxrange = [-352 -1]; %consider now the ima_k_spa_data is oversampled in kx
TSE.kyrange = [-176 -1]; %consider now the ima_k_spa_data is oversampled in kx

TSE.dyn_dim = dyn_nr;

%parameters for DPsti_TSE_phase_error_cor
pars.sense_map = 'external';  % external or ecalib 

pars.data_fn = data_fn;
pars.sense_ref = sense_ref_fn;
pars.coil_survey = coil_survey_fn;
pars.nav_phase_sm_kernel = 3;  %3 or 5, 1:no soomthing

pars.enabled_ch = [5 8 16];
pars.b0_shots = 127:168; %[] means first dynamic
pars.nonb0_shots = [1:42] + (1-1)*42;
pars.recon_x_locs = 1:160;

%paraemter for msDWIrecon called by DPsti_TSE_phase_error_cor
pars.msDWIrecon = initial_msDWIrecon_Pars;
pars.msDWIrecon.CG_SENSE_I.lamda=1;
pars.msDWIrecon.CG_SENSE_I.nit=5;
pars.msDWIrecon.CG_SENSE_I.tol = 1e-10;
pars.msDWIrecon.POCS.Wsize = [15 15];  %no point to be bigger than navigator area
pars.msDWIrecon.POCS.nit = 50;
pars.msDWIrecon.POCS.tol = 1e-10;
pars.msDWIrecon.POCS.lamda = 1;
pars.msDWIrecon.POCS.nufft = false;

pars.msDWIrecon.method='CG_SENSE_I'; %POCS_ICE CG_SENSE_I CG_SENSE_K LRT

%------------sense mask calc----------%
os = [1, 1, 1];
dim = [range(TSE.kxrange), range(TSE.kyrange), range(TSE.kzrange) ]+1;
sense_map_temp = get_sense_map_external(pars.sense_ref, pars.data_fn, pars.coil_survey, [dim(2) dim(2) dim(3)], os);
rs_command = sprintf('resize -c 0 %d', dim(1));
sense_map_temp = bart(rs_command, sense_map_temp);

TSE.sense_mask = abs(sense_map_temp(:,:,:,1 ))>0;
TSE_sense_map = sense_map_temp; %[]; %calc again using get_sense_map_external
clear sense_map_temp;    
%-------------------end---------------%

clear mr nav_im_recon_nufft nav_im_recon_nufft_1dyn nav_k_spa_data ima_kspa_sorted ima_default_recon
image_corrected = DPsti_TSE_phase_error_cor(ima_k_spa_data, TSE, TSE_sense_map, nav_data, pars);
% TODO make DPsti_TSE_phase_error_cor for POCS_ICE option

%% DTI data processing. ADC, FA map
b = 800;
g_all = [0 0 0;...
    1.00000,  0.00000,  1.00000; ...
    -0.41413,  0.94518,  0.96702; ...
    -1.00000,  0.00000,  1.00000; ...
    -0.41413, -0.94518,  0.96702; ...
    -0.96702, -0.94518,  0.41413; ...
    -0.96702,  0.94518,  0.41413; ...
    0.78398, -0.87792,  0.78398; ...
    0.78398,  0.87792,  0.78398];

clear MD FA eigvec

selected_volume = 1:9;
DTI_data = abs(image_corrected_all(:,:,:,selected_volume));
g = g_all(selected_volume,:);

[MD, FA, eigvec] = DTI_fitting(DTI_data, g, b);


mask = DTI_data(:,:,:,1)>10;
MD = bsxfun(@times,MD, mask );
FA = bsxfun(@times,FA, mask );
eigvec = bsxfun(@times,eigvec, mask );

figure(63); 
imshow(MD,[0 0.003]); colormap jet; colorbar; title('MD');
figure(64); 
imshow(FA,[0 1]);  colorbar; title('FA'); colormap hot;

eigvec = permute(eigvec,[1 2 3 5 4]);
figure(65);montage(permute(squeeze(eigvec(:,:,1,:,:)),[1 2 3 4]));  colorbar; title('eigenvector #1');


%% Linear motion correction
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
rigid_motion_correction = 0;
if(rigid_motion_correction)
    %% unwrap nav phase (2D) & fast rigid motion estimation
    close all; clc
    disp(' unwrap nav phase (2D) & fast rigid motion estimation...')
    
    clear PE_estimation
    
    ref_shot_ix = 1;
    valid_data_ratio = 0.4; % pixels with highest 40% intensity are used for processing
    
    nav_ima_phase_unwrapped = [];
    for dyn = 1:dyn_nr
        
        nav_im_to_unwrap = permute(nav_im_recon_nufft(:,:,1,:,:,dyn), [1 2 4 5 6 3]);  %[x,y,z,ch,shot,dyn] ->[x,y,shot,ch,dyn,z] this is 2D
        nav_ima_phase_unwrapped_dyn = spiral_nav_phase_unwrapping_2D(nav_im_to_unwrap, ref_shot_ix);
        nav_ima_phase_unwrapped = cat(4,nav_ima_phase_unwrapped, nav_ima_phase_unwrapped_dyn);
        
        [nav_ima_phase_unwrapped_diff, fitted_nav_ima_phase, linear_phase_xy,global_phase, global_phase_diff_initial] = ...
            rigidMotion_parameter_calculation(nav_im_to_unwrap, nav_ima_phase_unwrapped_dyn, ref_shot_ix, valid_data_ratio);  %nav_kspa_to_process  in kspace and with size of [kx ky ch shot]
        
        PE_estimation.linear_phase_xy_all(:,:,:,dyn) = (linear_phase_xy);
        PE_estimation.global_phase_all(:,:,dyn) = (global_phase);
        
    end
    PE_estimation.nav_kx_dim = size(nav_im_recon_nufft, 1);
    PE_estimation.nav_ky_dim = size(nav_im_recon_nufft, 2);
    
    disp('-finished- ');
    
    
    %% TSE data linear phase error correction
    correction_shot_range = 15:28;
    
    
    raw_fn.sense_ref_fn = sense_ref_fn;
    raw_fn.data_fn = data_fn;
    raw_fn.coil_survey_fn = coil_survey_fn;
    im_recon_nufft_cor = Perform_2D_SN_DPsti_recon(ima_k_spa_data, TSE, PE_estimation,  correction_shot_range, raw_fn);
end