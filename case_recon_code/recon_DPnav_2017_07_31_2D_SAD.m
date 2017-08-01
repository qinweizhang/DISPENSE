
clear;clc; close all;
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_07_31_SoSAD_SENSE_BRAIN/2017_07_31/DP_46368');
%% trajectory calculation
trj_save_fn('traj_Sc5_6_7.mat');
spira_3D_trjectory_calculation(trj_save_fn, 15);


%% TSE image and spiral Nav. data loading
clear
fn = 'dp_31072017_1954081_11_2_wip_dpsti_sosad_linearV4.raw';
save_fn = 'data_Sc11_3D.mat';

nav_kspa_data_read(fn, save_fn);
ima_kspa_data(fn, save_fn);


%% -----BART recon -------%
clear; close all; 
data_fn = 'data_Sc02_3D.mat';
trj_fn = 'traj_Sc5_6_7.mat';

%======recon parameters ======%
bart_recon.ignor_kz = 0;  % 1for yes. 0 for no
bart_recon.diffusion_nr = 1;
bart_recon.skip_point =0;
bart_recon.end_point = [];
bart_recon.recon_dim = [64, 64, 32];
bart_recon.trj_scale_dim = [24, 24, 15];
bart_recon.shot_nr = 1;
bart_recon.nsa_nr = 1;
bart_recon.sense_ref = 'dp_31072017_1915255_1000_7_wip_senserefscanV4.raw';
bart_recon.data_fn = 'dp_31072017_1915568_2_2_wip_dpsti_sosad_linearV4.raw';
bart_recon.coil_survey = 'dp_31072017_1911515_1000_2_wip_coilsurveyscanV4.raw';
bart_recon.update_SENSE_map = 0;
bart_recon.PICS = true;
bart_recon.sense_calc_method = 'external'; %'ecalib' or 'external'
%========end==================%
[reco_pics, igrid, igrid_rss] = bart_nufft_recon(data_fn, trj_fn, bart_recon);
 
figure(204); montage(permute(abs(reco_pics),[1 2 4 3]),'displayrange',[])
figure(205); montage(permute(angle(reco_pics),[1 2 4 3]),'displayrange',[-pi pi]); colormap jet
save(data_fn, 'reco_pics','igrid','igrid_rss','-append');

%% NUFFT recon. 
clear; close all; clc;

data_fn = 'data_Sc11_3D.mat';
trj_fn = 'traj_Sc5_6_7.mat'; 



%=============== recon parameters =========================
recon_par.ignore_kz = 1;
recon_par.recon_dim  = [24 24 1];
recon_par.dyn_nr = 1;
recon_par.skip_point = 0;
recon_par.end_point = 2000; %or []: till the end; 
recon_par.sense_map_recon = 0;
recon_par.update_SENSE_map = 0;
recon_par.sense_calc_method = 'external'; %'ecalib' or 'external'
recon_par.data_fn = 'dp_31072017_1954081_11_2_wip_dpsti_sosad_linearV4.raw';
recon_par.sense_ref = 'dp_31072017_1915255_1000_7_wip_senserefscanV4.raw';
recon_par.coil_survey = 'dp_31072017_1911515_1000_2_wip_coilsurveyscanV4.raw';
%========================  END  =========================

nav_im_recon_nufft = NUFFT_3D_recon(data_fn,trj_fn,recon_par);
save(data_fn, 'nav_im_recon_nufft','-append');
