
clear;clc; close all;
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_08_03_SoSAD_phantom');
%% trajectory calculation
close all; clear; clc;
trj_save_fn = 'traj_Sc16_17_18.mat';
trajectory_measure_distance = 15; %in mm
spira_3D_trjectory_calculation(trj_save_fn, trajectory_measure_distance);


%% TSE image and spiral Nav. data loading
clear
fn = 'dp_03082017_1849161_15_2_wipsc2dpstisosadlinearV4.raw';
save_fn = 'data_Sc15_3D.mat';

nav_kspa_data_read(fn, save_fn);
% ima_kspa_data(fn, save_fn);


%% NUFFT recon. 
clear; close all; clc;

data_fn = 'data_Sc15_3D.mat';
trj_fn = 'traj_Sc16_17_18.mat'; 

%=============== recon parameters =========================
recon_par.ignore_kz = 0;
recon_par.recon_dim  = [50 50 16];
recon_par.dyn_nr = 1;
recon_par.skip_point = 0;
recon_par.end_point = []; %or []: till the end; 
recon_par.interations = 10;
recon_par.sense_map_recon = 1;
recon_par.update_SENSE_map = 1;
recon_par.sense_calc_method = 'external'; %'ecalib' or 'external'
recon_par.data_fn = 'dp_03082017_1849161_15_2_wipsc2dpstisosadlinearV4.raw';
recon_par.sense_ref = 'dp_03082017_1810092_1000_7_wipsenserefscanV4.raw';
recon_par.coil_survey = 'dp_03082017_1806457_1000_2_wipcoilsurveyscanV4.raw';
%========================  END  =========================

nav_im_recon_nufft = NUFFT_3D_recon(data_fn,trj_fn,recon_par);
save(data_fn, 'nav_im_recon_nufft','-append');

%% -----BART recon -------%
clear; close all; 
data_fn = 'data_Sc05_3D.mat';
trj_fn = 'traj_Sc7_8_9.mat';

%======recon parameters ======%
bart_recon.ignor_kz = 0;  % 1for yes. 0 for no
bart_recon.diffusion_nr = 1;
bart_recon.skip_point =0;
bart_recon.end_point = [];
bart_recon.recon_dim = [50, 50, 33];
bart_recon.trj_scale_dim = [20, 20, 11];
bart_recon.shot_nr = 1;
bart_recon.nsa_nr = 1;
bart_recon.data_fn = 'dp_03082017_1835019_11_2_wipsc2dpstisosadlinearV4.raw';
bart_recon.sense_ref = 'dp_02082017_1516422_1000_7_wip_senserefscanV4.raw';
bart_recon.coil_survey = 'dp_02082017_1513421_1000_2_wip_coilsurveyscanV4.raw';
bart_recon.update_SENSE_map = 0;
bart_recon.PICS = true;
bart_recon.sense_calc_method = 'ecalib'; %'ecalib' or 'external'
%========end==================%
[reco_pics, igrid, igrid_rss] = bart_nufft_recon(data_fn, trj_fn, bart_recon);
 
figure(204); montage(permute(abs(reco_pics),[1 2 4 3]),'displayrange',[])
figure(205); montage(permute(angle(reco_pics),[1 2 4 3]),'displayrange',[-pi pi]); colormap jet
save(data_fn, 'reco_pics','igrid','igrid_rss','-append');

