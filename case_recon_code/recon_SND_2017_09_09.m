
clear; clc; close all
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_09_09_SND/watermelon')
%% trajectory calculation
close all; clear; clc;
trj_save_fn = 'traj_for_Sc20.mat';
trajectory_measure_distance = 15; %in mm
spira_3D_trjectory_calculation(trj_save_fn, trajectory_measure_distance);
disp('-finished- ');

%% SET path for all the following steps
clear; close all; clc

data_fn = '3d_09092017_1742514_25_2_wipsc17dpstisosadlinearV4.raw';
sense_ref_fn = '3d_09092017_1727115_1000_28_wipsenserefscanV4.raw';
coil_survey_fn  = '3d_09092017_1725270_1000_23_wipcoilsurveyscanV4.raw';

trj_mat_fn = 'traj_for_Sc19.mat';

%% Spiral Nav. data loading
disp('spiral Nav. data loading...')
nav_k_spa_data = nav_kspa_data_read(data_fn);

disp('-finished- ');
%% Spiral NUFFT recon.
disp(' Spiral NUFFT recon...');
close all;
[kx_length ch_nr shot_nr dyn_nr] = size(nav_k_spa_data);

nav_im_recon_nufft = [];
for dyn = 1:dyn_nr
    %=============== recon parameters =========================
    recon_par.ignore_kz = 0;
    recon_par.recon_dim  = [26 26 10];
    recon_par.dyn_nr = dyn;
    recon_par.skip_point = 0 ;
    recon_par.end_point = []; %or []: till the end;
    recon_par.interations = 10;
    recon_par.lamda = 0.1;
    recon_par.recon_all_shot = 1;
    recon_par.sense_map_recon = 1;
    recon_par.update_SENSE_map = 0;
    recon_par.sense_calc_method = 'external'; %'ecalib' or 'external'
    recon_par.data_fn = data_fn;
    recon_par.sense_ref = sense_ref_fn;
    recon_par.coil_survey = coil_survey_fn;
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
    
    
    nav_im_recon_nufft_1dyn = NUFFT_3D_recon(nav_k_spa_data,trj_mat_fn,recon_par, nav_sense_map);
    nav_im_recon_nufft = cat(6, nav_im_recon_nufft, nav_im_recon_nufft_1dyn);
end
% nav_sense_map = circshift(nav_sense_map, round(17.26/115.00*size(nav_sense_map,1)));
% nav_im_recon_nufft = circshift(nav_im_recon_nufft, -1*round(17.26/115.00*size(nav_sense_map,1)));

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

[ima_k_spa_data,TSE.ky_matched,TSE.kz_matched,TSE.shot_matched, TSE.ch_dim,ima_kspa_sorted, ima_default_recon, TSE_sense_map, TSE.kxrange, TSE.kyrange, TSE.kzrange] = ...
    TSE_data_sortting(data_fn, sense_ref_fn, coil_survey_fn);


figure(606); immontage4D(permute(abs(ima_default_recon(80:240,:,:,:)),[1 2 4 3]), [10 120]);
disp('-finished- ');

%% TSE data non-rigid phase error correction (iterative) CG_SENSE
nav_data = reshape(nav_im_recon_nufft, size(nav_im_recon_nufft,1), size(nav_im_recon_nufft, 2), size(nav_im_recon_nufft, 3), max(TSE.shot_matched));

TSE.kxrange = [-320 -1];
TSE.dyn_dim = dyn_nr;
TSE_sense_map = []; %calc again using get_sense_map_external

%parameters for DPsti_TSE_phase_error_cor
pars.sense_map = 'external';  % external or ecalib

pars.data_fn = data_fn;
pars.sense_ref = sense_ref_fn;
pars.coil_survey = coil_survey_fn;

pars.enabled_ch = [1:2];
pars.b0_shots = []; %[] means first dynamic
pars.nonb0_shots = 28:54;
pars.recon_x_locs = 120:220;

%paraemter for msDWIrecon called by DPsti_TSE_phase_error_cor
pars.msDWIrecon = initial_msDWIrecon_Pars;
pars.msDWIrecon.CG_SENSE_I.lamda=0.1;
pars.msDWIrecon.CG_SENSE_I.nit=10;
pars.msDWIrecon.CG_SENSE_I.tol = 1e-10;
pars.msDWIrecon.POCS.Wsize = [15 15];  %no point to be bigger than navigator area
pars.msDWIrecon.POCS.nit = 50;
pars.msDWIrecon.POCS.tol = 1e-10;
pars.msDWIrecon.POCS.lamda = 1;
pars.msDWIrecon.POCS.nufft = false;

pars.msDWIrecon.method='POCS_ICE'; %POCS_ICE CG_SENSE_I CG_SENSE_K LRT

clear mr nav_sense_map nav_im_recon_nufft nav_im_recon_nufft_1dyn nav_k_spa_data ima_kspa_sorted ima_default_recon
image_corrected = DPsti_TSE_phase_error_cor(ima_k_spa_data, TSE, TSE_sense_map, nav_data, pars);
% TODO make DPsti_TSE_phase_error_cor for POCS_ICE option


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