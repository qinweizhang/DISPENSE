
clear; clc; close all
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2018_01_28_SND_brain')
%% trajectory calculation
close all; clear; clc;
trj_save_fn = 'traj_for_Sc2.mat';
trajectory_measure_distance = 15; %in mm
spira_3D_trjectory_calculation(trj_save_fn, trajectory_measure_distance);
disp('-finished- ');

%% SET path for all the following steps
clear; close all; clc

data_fn = 'sn_28012018_1246416_2_2_wip_sc2_3d_snd_brain_4bV4.raw';
sense_ref_fn = 'sn_28012018_1246096_1000_5_wip_senserefscanV4.raw';
coil_survey_fn  = 'sn_28012018_1243593_1000_2_wip_coilsurveyscanV4.raw';

trj_mat_fn = 'traj_for_Sc2_old.mat';

%% Spiral Nav. data loading
disp('spiral Nav. data loading...')
coil_cmp_nr = 20;
[nav_k_spa_data, Nav_VirtualCoilMartix] = nav_kspa_data_read(data_fn, coil_cmp_nr);
% nav_k_spa_data = nav_kspa_data_read(data_fn);

disp('-finished- ');
%% Spiral NUFFT recon.
disp(' Spiral NUFFT recon...');
save_mat_fn = 'Sc02.mat';
close all;
[kx_length ch_nr shot_nr, dyn_nr] = size(nav_k_spa_data);

offcenter_xy = [0 0]; 
FOV_xy = [250 159.0909 ];
% nav_im_recon_nufft = [];
dyn_recon = 9:-1:1;
for d = 1:length(dyn_recon)
    tic
    dyn  = dyn_recon(d);
    disp(['dynamic: ',num2str(dyn)]);
    %=============== recon parameters =========================
    recon_par.ignore_kz = 0;
    recon_par.acq_dim = [32 32 18];  
    recon_par.recon_dim  = [32 32 18];
    recon_par.dyn_nr = dyn;
    recon_par.skip_point = 0 ;
    recon_par.end_point = []; %or []: till the end;
    recon_par.interations = 10;
    recon_par.lamda = 0;
    recon_par.recon_all_shot = 1;
    recon_par.sense_map_recon =1; 
    recon_par.update_SENSE_map = 0;
    recon_par.sense_calc_method = 'external'; %'ecalib' or 'external'
    recon_par.sense_os = [1 FOV_xy(1)/FOV_xy(2)];  %oversampling in x and y: control sense FOV
    recon_par.data_fn = data_fn;
    recon_par.sense_ref = sense_ref_fn;
    recon_par.coil_survey = coil_survey_fn;
    recon_par.parfor = 1; %parfor recon for different shots
    
    recon_par.channel_by_channel = 1;
    recon_par.channel_by_channel = recon_par.channel_by_channel .* (1-recon_par.sense_map_recon );
    %========================  END  =========================
     if(~exist('nav_sense_map', 'var')&&recon_par.sense_map_recon)
        recon_par.update_SENSE_map = 1;
     end
    
    if(recon_par.update_SENSE_map)
        clc_sens_par.cc = 1;
        clc_sens_par.disp = 0;
        [nav_sense_map, nav_sense_Psi] = calc_sense_map(recon_par.data_fn, recon_par.sense_ref,  recon_par.coil_survey, recon_par.recon_dim,recon_par.sense_calc_method, recon_par.sense_os, clc_sens_par);
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
    save(save_mat_fn, 'nav_im_recon_nufft', 'dyn_nr', '-append'); 
    toc
end
% nav_sense_map = circshift(nav_sense_map, round(17.26/115.00*size(nav_sense_map,1)));
% nav_im_recon_nufft = circshift(nav_im_recon_nufft, -1*round(17.26/115.00*size(nav_sense_map,1)));
dyn = 2;
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

%% TSE data sorting and default recon
close all; clc; clear TSE
disp(' TSE data sorting and default recon...')

parameter2read.dyn = [];
parameter2read.cc_nr = 20; %20;
parameter2read.sense_recon = 0;

[ima_k_spa_data,TSE.ky_matched,TSE.kz_matched,TSE.shot_matched, TSE.ch_dim,ima_kspa_sorted, ima_default_recon, TSE_sense_map, TSE.kxrange, TSE.kyrange, TSE.kzrange, TSE.VirtualCoilMartix] = ...
    TSE_data_sortting(data_fn, sense_ref_fn, coil_survey_fn,parameter2read);

figure(610); immontage4D(permute(abs(ima_default_recon(:,:,:,:)),[1 2 4 3]), []);

TSE
assert(length(TSE.ky_matched)==size(ima_k_spa_data,2),'Profile number does not match with data size!')
disp('-finished- ');

%% TSE data non-rigid phase error correction (iterative) CG_SENSE
save_mat_fn = 'Sc02.mat';

nav_data = reshape(nav_im_recon_nufft, size(nav_im_recon_nufft,1), size(nav_im_recon_nufft, 2), size(nav_im_recon_nufft, 3), max(TSE.shot_matched));

TSE.dyn_dim = dyn_nr;
TSE.SENSE_kx =1;
TSE.SENSE_ky =2;
TSE.SENSE_kz =1;

% TSE.kyrange = [-38 -1]; 
TSE.kxrange = [-352 -1]; %consider now the ima_k_spa_data is oversampled in kx; kx oversmapled by 2 + 
TSE.kzrange = [-68, -1];

TSE.Ixrange = [ceil(TSE.kxrange(1).*TSE.SENSE_kx) -1];
TSE.Iyrange = [ceil(TSE.kyrange(1).*TSE.SENSE_ky) -1];
TSE.Izrange = [ceil(TSE.kzrange(1).*TSE.SENSE_kz) -1];
TSE.kyrange = TSE.Iyrange;
TSE.kzrange = TSE.Izrange;



%parameters for DPsti_TSE_phase_error_cor
pars.sense_map = 'external';  % external or ecalib 

pars.data_fn = data_fn;
pars.sense_ref = sense_ref_fn;
pars.coil_survey = coil_survey_fn;
pars.parfor = 1;
pars.nav_phase_sm_kernel = 3;  %3 or 5, 1:no soomthing
pars.recon_x_locs = 88:(88+176);
pars.enabled_ch = 1:TSE.ch_dim;
pars.b0_shots = []; %[] means first dynamic
pars.recon_dyn = 9; %9:-1:1;
pars.large_scale_recon = 0; %true; % Choose to use DPsti_TSE_phase_error_cor_large_scale.m or DPsti_TSE_phase_error_cor.m
pars.nocorrection = 0;


%paraemter for msDWIrecon called by DPsti_TSE_phase_error_cor
pars.msDWIrecon = initial_msDWIrecon_Pars;
pars.msDWIrecon.trim_kspa_filter_mask_size = 7;
pars.msDWIrecon.CG_SENSE_I.lamda=1e-2;
pars.msDWIrecon.CG_SENSE_I.nit=15;
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

clc_sens_par.cc = 1;
clc_sens_par.disp = 0;
[sense_map_temp, TSE.sense_Psi] = get_sense_map_external(pars.sense_ref, pars.data_fn, pars.coil_survey, [dim(1)/2 dim(2) dim(3)], os, clc_sens_par);
%----compress sense map and sense_Psi
if(isfield(TSE, 'VirtualCoilMartix'))
    if(~isempty(TSE.VirtualCoilMartix))
        [sense_map_temp, TSE.sense_Psi] = compress_sense_map_Psi(TSE.VirtualCoilMartix, sense_map_temp,  TSE.sense_Psi);
    end
end
        
rs_command = sprintf('resize -c 0 %d', dim(1));
sense_map_temp = bart(rs_command, sense_map_temp);

TSE.sense_mask = abs(sense_map_temp(:,:,:,1 ))>0;
TSE_sense_map = sense_map_temp; %[]; %calc again using get_sense_map_external
clear sense_map_temp;    
%-------------------end---------------%

clear mr nav_im_recon_nufft nav_im_recon_nufft_1dyn nav_k_spa_data ima_kspa_sorted ima_default_recon



shot_per_dyn = max(TSE.shot_matched) / TSE.dyn_dim;

for dx = 1:length(pars.recon_dyn)
    tic
    d = pars.recon_dyn(dx)
    
    pars.nonb0_shots = [1:shot_per_dyn] + (d-1)*shot_per_dyn;
    
    if(pars.large_scale_recon)
        if(~pars.nocorrection)
            result = DPsti_TSE_phase_error_cor_large_scale(ima_k_spa_data, TSE, TSE_sense_map, nav_data, pars);
        else
            result = DPsti_TSE_phase_error_cor_large_scale(ima_k_spa_data, TSE, TSE_sense_map, ones(size(nav_data)), pars);
        end
    else
        if(~pars.nocorrection)
            result = DPsti_TSE_phase_error_cor(ima_k_spa_data, TSE, TSE_sense_map, nav_data, pars);
        else
            result = DPsti_TSE_phase_error_cor(ima_k_spa_data, TSE, TSE_sense_map, ones(size(nav_data)), pars);
        end
    end

    if(~pars.nocorrection)
        image_corrected(:,:,:,d)  = result;
        save(save_mat_fn,'image_corrected','-append');
    else
        image_uncorrected(:,:,:,d)  = result;
        save(save_mat_fn,'image_uncorrected','-append');
    end
     
     
     
     elaps_t=toc;
     if(~pars.parfor)
         msg = sprintf(['Recon finishted for {', data_fn,'} ; dynamic %d ; duration %f; s', 10, 'Saved as ', save_mat_fn],d, elaps_t);
         sendmail_from_yahoo('q.zhang@amc.nl','Matlab Message',msg);
     end
end
% TODO make DPsti_TSE_phase_error_cor for POCS_ICE option

%% DTI data processing. ADC, FA map
% b = 800;
% 
% 
% g_all = [0.000,  0.000, -0.668,  0.000,  0.668,  0.668, -0.668,  0.621, -0.621;...
%          0.000,  0.707, -0.293, -0.707, -0.293, -0.684, -0.684,  0.554,  0.554;...
%          0.000,  0.707,  0.684,  0.707,  0.684,  0.293,  0.293,  0.554,  0.554];
% g_all = g_all';
% 
% clear MD FA eigvec
% 
% selected_volume = 1:9;
% slice = 20;
% DTI_data = abs(image_corrected(:,:,slice,selected_volume));
% g = g_all(selected_volume,:);
% 
% [MD, FA, eigvec] = DTI_fitting(DTI_data, g, b);
% 
% 
% mask = DTI_data(:,:,:,1)>10;
% MD = bsxfun(@times,MD, mask );
% FA = bsxfun(@times,FA, mask );
% eigvec = bsxfun(@times,eigvec, mask );
% 
% figure(63); 
% imshow(MD,[0 0.003]); colormap jet; colorbar; title('MD');
% figure(64); 
% imshow(FA,[0.2 1]);  colorbar; title('FA'); colormap hot;
% 
% eigvec = permute(eigvec,[1 2 3 5 4]);
% figure(65);montage(permute(squeeze(eigvec(:,:,1,:,:)),[1 2 3 4]));  colorbar; title('eigenvector #1');


