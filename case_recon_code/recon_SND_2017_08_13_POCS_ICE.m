clear;clc; close all;
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_08_13_SoSAD_phantom_pormelo');
%% SET path for all the following steps
clear; close all; clc

data_fn = 'sa_13082017_1958269_15_2_wip_sc18_dpsti_sosad_linearV4.raw';
sense_ref_fn = 'sa_13082017_1957331_1000_34_wip_senserefscanV4.raw';
coil_survey_fn  = 'sa_13082017_1956371_1000_29_wip_coilsurveyscanV4.raw';

trj_mat_fn = 'traj2_Sc29_27_28_for_Sc4.mat';


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
    recon_par.ignore_kz = 1;
    recon_par.recon_dim  = [36 36 1];
    recon_par.dyn_nr = dyn;
    recon_par.skip_point = 0 ;
    recon_par.end_point = 2000; %[]; %or []: till the end;
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
%% TSE data loading and default recon
close all; clc
disp(' TSE data sorting and default recon...')

[ima_k_spa_data,TSE.ky_matched,TSE.kz_matched,TSE.shot_matched, TSE.ch_dim,ima_kspa_sorted,...
    ima_default_recon, TSE_sense_map, TSE.kxrange, TSE.kyrange, TSE.kzrange] = ...
    TSE_data_sortting(data_fn, sense_ref_fn, coil_survey_fn);


figure(606); immontage4D(permute(abs(ima_default_recon(80:240,:,:,:)),[1 2 4 3]), [10 120]);
disp('-finished- ');

%% POCS-ICE

TSE.kxrange = [-320 -1];
TSE.dyn_dim = 4;
TSE_sense_map = []; %calc again using get_sense_map_external

nav_im = reshape(nav_im_recon_nufft, size(nav_im_recon_nufft,1), size(nav_im_recon_nufft, 2), size(nav_im_recon_nufft, 3), max(TSE.shot_matched));

%parameters for DPsti_TSE_phase_error_cor
pars.sense_map = 'external';  % external or ecalib

pars.data_fn = data_fn;
pars.sense_ref = sense_ref_fn;
pars.coil_survey = coil_survey_fn;
pars.nav_trj = trj_mat_fn;

pars.enabled_ch = [1:38];
pars.b0_shots = [];
pars.nonb0_shots = 15:56;
pars.recon_x_locs = 120:220;


%paraemter for msDWIrecon called by DPsti_TSE_phase_error_cor
pars.msDWIrecon = initial_msDWIrecon_Pars;
pars.msDWIrecon.CG_SENSE_I.lamda=0.1;
pars.msDWIrecon.CG_SENSE_I.nit=10;
pars.msDWIrecon.CG_SENSE_I.tol = 1e-10;
pars.msDWIrecon.method='POCS_ICE';

image_corrected = DPsti_TSE_phase_error_cor(ima_k_spa_data, TSE, TSE_sense_map, nav_k_spa_data, pars);





%%
tic;

k0_idx = [floor(kx_dim/2)+1 floor(ky_dim/2)+1 floor(kz_dim/2)+1];
for sh =1:sh_dim
    if(abs(kspa_xyz(k0_idx(1),k0_idx(2),k0_idx(3), 1, sh))>0)
        ref_shot = sh;
    end
end
phase_error_3D = bsxfun(@rdivide, phase_error_3D, phase_error_3D(:,:,:,ref_shot)); %difference with the first
phase_error_3D = permute(normalize_sense_map(squeeze(phase_error_3D)),[1 2 3 5 4]); %miss use normalize_sense_map

%=========select data. fixed=============================================
kspa = squeeze(kspa_x_yz(recon_x_loc, :, :, :, :));
sense_map = squeeze(sense_map_3D(recon_x_loc,:,:,:));
phase_error = permute(squeeze(phase_error_3D(recon_x_loc,:,:,:,:)),[1 2 4 3]);
%========================================================================

recon_x_loc = 90;

%recon parameter
pars.POCS.Wsize = [15 15];  %no point to be bigger than navigator area
pars.POCS.nit = 50;
pars.POCS.tol = 1e-10;
pars.POCS.lamda = 1;
pars.POCS.nufft = false;

pars.method='POCS_ICE'; %POCS_ICE CG_SENSE_I CG_SENSE_K LRT

image_corrected = msDWIrecon(kspa, (sense_map), (phase_error), pars);


%display
figure(102);
subplot(131);imshow(squeeze(abs(ima_ref_rss(recon_x_loc,:,:))),[]); title('reference');
subplot(132);imshow(squeeze(abs(im_recon_direct(recon_x_loc,:,:))),[]); title('direct recon');
subplot(133);imshow(squeeze(abs(image_corrected)),[]); title('msDWIrecon'); xlabel(pars.method);

toc;