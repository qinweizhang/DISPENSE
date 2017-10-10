
clear;clc; close all;
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_08_04_SoSAD_phantom');
%% trajectory calculation
close all; clear; clc;
trj_save_fn = 'traj_Sc16_17_18.mat';
trajectory_measure_distance = 15; %in mm
spira_3D_trjectory_calculation(trj_save_fn, trajectory_measure_distance);
disp('-finished- ');

%% SET path for all the following steps
clear; close all; clc

data_fn = 'dp_04082017_1341243_2_2_wip_sc6_dpsti_sosad_linear-ppuV4.raw';
sense_ref_fn = 'dp_04082017_1340473_1000_7_wip_senserefscanV4.raw';
coil_survey_fn  = 'dp_04082017_1336408_1000_2_wip_coilsurveyscanV4.raw';

trj_mat_fn = 'traj_Sc9_10_11_for_Sc2.mat';

%% Spiral Nav. data loading
disp('spiral Nav. data loading...')
nav_k_spa_data = nav_kspa_data_read(data_fn);

disp('-finished- ');
%% Spiral NUFFT recon.
disp(' Spiral NUFFT recon...');
close all;
[kx_length, ch_nr shot_nr, dyn_nr] = size(nav_k_spa_data);

nav_im_recon_nufft = [];
for dyn = 1:dyn_nr
    dyn
    %=============== recon parameters =========================
    recon_par.ignore_kz = 0;
    recon_par.acq_dim = [36 36 13];  
    recon_par.recon_dim  = [36 36 13];
    recon_par.dyn_nr = dyn;
    recon_par.skip_point = 0 ;
    recon_par.end_point = []; %or []: till the end;
    recon_par.interations = 10;
    recon_par.lamda = 0;
    recon_par.recon_all_shot = 0;
    recon_par.sense_map_recon =0; 
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
    
    nav_sense_Psi = [];
    if(recon_par.update_SENSE_map)
        [nav_sense_map, nav_sense_Psi] = calc_sense_map(recon_par.data_fn, recon_par.sense_ref,  recon_par.coil_survey, recon_par.recon_dim,recon_par.sense_calc_method);
    end
    
    if(recon_par.sense_map_recon == 0)
        nav_sense_map = ones([recon_par.recon_dim ch_nr]);
    end
    nav_sense_map = normalize_sense_map(nav_sense_map);
    
    offcenter_xy = [-4 6]; 
    FOV_xy = [220 100];
    
    nav_im_recon_nufft_1dyn = NUFFT_3D_recon(nav_k_spa_data,trj_mat_fn,recon_par, nav_sense_map, nav_sense_Psi,offcenter_xy, FOV_xy);
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

%% TSE data sorting and default recon
close all; clc
disp(' TSE data sorting and default recon...')

parameter2read.dyn = [];

[ima_k_spa_data,TSE.ky_matched,TSE.kz_matched,TSE.shot_matched, TSE.ch_dim,ima_kspa_sorted, ima_default_recon, TSE_sense_map, TSE.kxrange, TSE.kyrange, TSE.kzrange] = ...
    TSE_data_sortting(data_fn, sense_ref_fn, coil_survey_fn,parameter2read);

TSE

figure(610); immontage4D(permute(abs(ima_default_recon(:,:,:,:)),[1 2 4 3]), []);
disp('-finished- ');

%% TSE data non-rigid phase error correction (iterative) CG_SENSE
nav_data = reshape(nav_im_recon_nufft, size(nav_im_recon_nufft,1), size(nav_im_recon_nufft, 2), size(nav_im_recon_nufft, 3), max(TSE.shot_matched));

TSE.SENSE_kx =1;
TSE.SENSE_ky =2;
TSE.SENSE_kz =2;

TSE.kyrange = [-96 -1]; 
TSE.kxrange = [-384 -1]; %consider now the ima_k_spa_data is oversampled in kx; kx oversmapled by 2 + 

TSE.Ixrange = [ceil(TSE.kxrange(1).*TSE.SENSE_kx) -1];
TSE.Iyrange = [ceil(TSE.kyrange(1).*TSE.SENSE_ky) -1];
TSE.Izrange = [ceil(TSE.kzrange(1).*TSE.SENSE_kz) -1];
TSE.kyrange = TSE.Iyrange;
TSE.kzrange = TSE.Izrange;

TSE.dyn_dim = dyn_nr;

%parameters for DPsti_TSE_phase_error_cor
pars.sense_map = 'external';  % external or ecalib 

pars.data_fn = data_fn;
pars.sense_ref = sense_ref_fn;
pars.coil_survey = coil_survey_fn;
pars.nav_phase_sm_kernel = 3;  %3 or 5, 1:no soomthing
pars.recon_x_locs = 96:258;
pars.enabled_ch = [7: 13];
pars.b0_shots = []; %[] means first dynamic


%paraemter for msDWIrecon called by DPsti_TSE_phase_error_cor
pars.msDWIrecon = initial_msDWIrecon_Pars;
pars.msDWIrecon.CG_SENSE_I.lamda=1e-2;
pars.msDWIrecon.CG_SENSE_I.nit=40;
pars.msDWIrecon.CG_SENSE_I.tol = 1e-10;
pars.msDWIrecon.POCS.Wsize = [15 15];  %no point to be bigger than navigator area
pars.msDWIrecon.POCS.nit = 50;
pars.msDWIrecon.POCS.tol = 1e-10;
pars.msDWIrecon.POCS.lamda = 1;
pars.msDWIrecon.POCS.nufft = false;

pars.msDWIrecon.method='CG_SENSE_I'; %POCS_ICE CG_SENSE_I CG_SENSE_K LRT

%------------sense mask calc----------%
os = [1, 1, 1];
dim = [range(TSE.Ixrange), range(TSE.Iyrange), range(TSE.Izrange) ]+1;
[sense_map_temp, TSE.sense_Psi] = get_sense_map_external(pars.sense_ref, pars.data_fn, pars.coil_survey, [dim(2) dim(2) dim(3)], os);
rs_command = sprintf('resize -c 0 %d', dim(1));
sense_map_temp = bart(rs_command, sense_map_temp);

TSE.sense_mask = abs(sense_map_temp(:,:,:,1 ))>0;
TSE_sense_map = sense_map_temp; %[]; %calc again using get_sense_map_external
clear sense_map_temp;    
%-------------------end---------------%

clear mr nav_im_recon_nufft nav_im_recon_nufft_1dyn nav_k_spa_data ima_kspa_sorted ima_default_recon

for d = 6:6
    pars.nonb0_shots = [1:40] + (d-1)*40;
    image_corrected(:,:,:,d) = DPsti_TSE_phase_error_cor(ima_k_spa_data, TSE, TSE_sense_map, nav_data, pars);
    save('Sc4.mat','image_sense_corrected','-append');
end
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


