
clear; clc; close all
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_11_05_SND')
%% trajectory calculation (OPTIONAL)
close all; clear; clc;
trj_save_fn = 'traj_for_Sc8.mat';
trajectory_measure_distance = 15; %in mm
spira_3D_trjectory_calculation(trj_save_fn, trajectory_measure_distance);
disp('-finished- ');

%% SET path for all the following steps
clear; close all; clc

data_fn = 'sn_05112017_1226064_2_2_wip_2d_snd_brain_3b_lrt_csV4.raw';
sense_ref_fn = 'sn_05112017_1225311_1000_5_wip_senserefscanV4.raw';
coil_survey_fn  = 'sn_05112017_1223074_1000_2_wip_coilsurveyscanV4.raw';

trj_mat_fn = 'traj_for_Sc2_5.mat';

%% Spiral Nav. data loading
disp('spiral Nav. data loading...')
[nav_k_spa_data, Nav_VirtualCoilMartix] = nav_kspa_data_read(data_fn);
% nav_k_spa_data = nav_kspa_data_read(data_fn);

disp('-finished- ');



%% TSE imaging data loading (default recon)
close all; clc
disp(' TSE data sorting and default recon...')

parameter2read.dyn = [];

[ima_k_spa_data,TSE.ky_matched,TSE.kz_matched,TSE.shot_matched, TSE.ch_dim,ima_kspa_sorted, ima_default_recon, TSE_sense_map, TSE.kxrange, TSE.kyrange, TSE.kzrange, TSE.VirtualCoilMartix] = ...
    TSE_data_sortting(data_fn, sense_ref_fn, coil_survey_fn,parameter2read);

figure(610); immontage4D(permute(abs(ima_default_recon(:,:,:,:)),[1 2 4 3]), [0 80]);

TSE
assert(length(TSE.ky_matched)==size(ima_k_spa_data,2),'Profile number does not match with data size!')
disp('-finished- ');


%% SET parameter

%-------------TSE pars.-------------%
TSE.kxrange = [-512 -1];
TSE.dyn_dim = size(ima_kspa_sorted, 4);
TSE.shot_per_dyn = max(TSE.shot_matched) / TSE.dyn_dim;

TSE.kx_dim = TSE.kxrange(2) - TSE.kxrange(1) + 1;
TSE.ky_dim = TSE.kyrange(2) - TSE.kyrange(1) + 1;  %max_ky * 2 + 1;
TSE.kz_dim = TSE.kzrange(2) - TSE.kzrange(1) + 1;
TSE.SENSE_kx =1;
TSE.SENSE_ky =1;
TSE.SENSE_kz =1;
TSE.Ixrange = TSE.kxrange;
TSE.Iyrange = TSE.kyrange;
TSE.Izrange = TSE.kzrange;

%------------- General pars.-------------%

save_mat_fn = 'Sc02.mat';

pars = initial_msDWIrecon_Pars;
pars.sense_map = 'external';  % external or ecalib 
pars.data_fn = data_fn;
pars.sense_ref = sense_ref_fn;
pars.coil_survey = coil_survey_fn;
pars.nav_phase_sm_kernel = 3;  %3 or 5, 1:no soomthing
pars.recon_x_locs = 120:400; %80:270;
pars.enabled_ch = 1:TSE.ch_dim;
pars.b0_shots = []; %[] means first dynamic
pars.nonb0_shots = 26:50;
if(isempty(pars.b0_shots))
    pars.b0_shots = 1:TSE.shot_per_dyn;
end            


%% Get TSE SENSE maps

%------------sense mask calc----------%
os = [1, 1, 1];
dim = [range(TSE.Ixrange), range(TSE.Iyrange), range(TSE.Izrange) ]+1;
[tse_sense_map_unpaded, TSE.sense_Psi] = get_sense_map_external(pars.sense_ref, pars.data_fn, pars.coil_survey, [dim(1)/2 dim(2) dim(3)], os);
%----compress sense map and sense_Psi
if(isfield(TSE, 'VirtualCoilMartix'))
    if(~isempty(TSE.VirtualCoilMartix))
        [tse_sense_map_unpaded, TSE.sense_Psi] = compress_sense_map_Psi(TSE.VirtualCoilMartix, tse_sense_map_unpaded,  TSE.sense_Psi);
    end
end
        
rs_command = sprintf('resize -c 0 %d', dim(1));
tse_sense_map_paded = bart(rs_command, tse_sense_map_unpaded);

TSE.sense_mask = abs(tse_sense_map_paded(:,:,:,1 ))>0;
TSE_sense_map = tse_sense_map_paded; %[]; %calc again using get_sense_map_external
figure(3); montage(abs(TSE_sense_map),'displayrange',[]); xlabel('SENSE maps')
%-------------------end---------------%

%% Get b0 data

%=======TSE imaging=======%
b0_kpa = sort_k_spa_sh_by_sh_beta(ima_k_spa_data, pars.b0_shots, TSE, pars);
%----combine shots directly; non-zero average in the shot dim
b0_kspa = sum(b0_kpa, 5)./ sum(abs(b0_kpa)>0, 5); b0_kspa(find(isnan(b0_kspa)+isinf(b0_kspa))) = 0;
b0_ima = bart('pics -l1 -r 0.5 ', b0_kspa, TSE_sense_map);
figure(2); imshow(abs(b0_ima),[]); xlabel('CS recon of b0 data: shot directly combined')
b0_kspa_full = squeeze(fft2d(bsxfun(@times, b0_ima, TSE_sense_map)));

%=======Nav=======%
nav_b0_kspa = mean(nav_k_spa_data(:,:,pars.b0_shots),3);  %average of all shots
% nav_b0_kspa = nav_k_spa_data(:,:,pars.b0_shots(end));  %use one shot

size(b0_kspa_full) 
size(nav_b0_kspa)
clear b0_kpa b0_kspa b0_ima



%% Get non-b0 data

%=======TSE imaging=======%
nonb0_kpa = squeeze(sort_k_spa_sh_by_sh_beta(ima_k_spa_data, pars.nonb0_shots, TSE, pars));

%=======Nav=======%
nav_nonb0_kspa = squeeze(nav_k_spa_data(:,:,pars.nonb0_shots));

size(nonb0_kpa) 
size(nav_nonb0_kspa)



%% Get NAV SENSE maps  + SET more parameters 

%------------sense mask calc----------%
nav_dim = [42 42 1];
FOV = [250 250];
nav_sense_os = [1 FOV_xy(1)/FOV_xy(2)];
[nav_sense_map, TSE.sense_Psi] = get_sense_map_external(pars.sense_ref, pars.data_fn, pars.coil_survey, [nav_dim(1) nav_dim(2) nav_dim(3)], nav_sense_os);
%----compress sense map and sense_Psi
if(isfield(TSE, 'VirtualCoilMartix'))
    if(~isempty(TSE.VirtualCoilMartix))
        [nav_sense_map, TSE.sense_Psi] = compress_sense_map_Psi(TSE.VirtualCoilMartix, nav_sense_map,  TSE.sense_Psi);
    end
end

figure(4); montage(abs(nav_sense_map),'displayrange',[]); xlabel('NAV SENSE maps')
%-------------------end---------------%


%----------------NUFFT pars. ----------------%
load(trj_mat_fn);
chosen_sense_map = nav_sense_map;  % or tse_sense_map_unpaded nav_sense_map
recon_par.recon_dim = [size(chosen_sense_map,1) size(chosen_sense_map,2)];
recon_par.acq_dim = [42 42];  %<<<<<<<<<<<<input this
recon_par.rescale_factor = recon_par.recon_dim./ recon_par.acq_dim;

scale_foctor_xy = max(2*pi/(max(trj_meas_kx)-min(trj_meas_kx))./recon_par.rescale_factor(1),2*pi/(max(trj_meas_ky)-min(trj_meas_ky))./recon_par.rescale_factor(2));
trj_nufft = [trj_meas_kx trj_meas_ky] * scale_foctor_xy;

pars.LRT.mix_trajectory = 1;
pars.LRT.NUFFT_nav_sense = nuFTOperator(trj_nufft,[size(chosen_sense_map,1) size(chosen_sense_map,2)],squeeze(chosen_sense_map),6);
pars.LRT.NUFFT_nav_1ch = nuFTOperator(trj_nufft,[size(chosen_sense_map,1) size(chosen_sense_map,2)],ones(size(chosen_sense_map,1), size(chosen_sense_map,2)),6);
pars.LRT.trj_length = length(trj_nufft);
pars.LRT.nav_dim = recon_par.acq_dim;

%-------------try one nufft recon------------%
tic
image = regularizedReconstruction(pars.LRT.NUFFT_nav_sense,col(nav_b0_kspa(:,:,1)),@L2Norm,0.5,'maxit',10,'tol', 1e-10);   %nav_b0_kspa nav_nonb0_kspa
k_temp = fft2d(image);

upsampling = 1; downsampling = ~upsampling;

if(upsampling)
    rs_command = sprintf('resize -c 0 %d 1 %d', size(tse_sense_map_unpaded,1), size(tse_sense_map_unpaded,2));
    k_interpo = bart(rs_command, k_temp);
    image_interpo = ifft2d(k_interpo).*(abs(tse_sense_map_unpaded(:,:,1))>0);
    
elseif(downsampling)  
    rs_command = sprintf('resize -c 0 %d 1 %d', size(nav_sense_map,1), size(nav_sense_map,2));
    k_interpo = bart(rs_command, k_temp);
    image_interpo = ifft2d(k_interpo).*(abs(nav_sense_map(:,:,1))>0);
end

figure(5); imagesc(abs(image_interpo)); axis off; axis equal; colormap gray
figure(6); imagesc(angle(image_interpo),[-pi pi]); axis off; axis equal; colormap jet

toc

%% RECON

clear kspa
%-----range all kspace data---%
for idx1 =1:length(pars.nonb0_shots)
    for idx2=1:2
        if(idx2==1) %nav colume
            kspa{idx1,idx2} =  nav_nonb0_kspa(:,:,idx1);
        else %tse colume
            kspa{idx1,idx2} =  nonb0_kpa(:,:,:,idx1);
        end
    end
end
kspa_b0{1,1} =  nav_b0_kspa;
kspa_b0{1,2} =  b0_kspa_full;
kspa = cat(1, kspa_b0, kspa);


%------LRT pars---------%
pars.method='LRT'; %POCS_ICE CG_SENSE_I CG_SENSE_K LRT

pars.LRT.Lg=3;
pars.LRT.L3=4;
pars.LRT.L4=1;
pars.LRT.mu = 2e2;
pars.LRT.beta = 1;
pars.LRT.lambda = 4e-2;

pars.LRT.sparsity_transform='TV';
%     pars.LRT.Imref=cat(3, repmat(squeeze(im_b0_ref(recon_x_loc,:,:,:)), [1 1 1 2]), repmat(squeeze(im_ref(recon_x_loc,:,:,:)), [1 1 length(recon_shot_range) 2]));
pars.LRT.x=88;  %loc in dim2
pars.LRT.y=300; %loc in dim1
pars.LRT.increase_penalty_parameters=false;
pars.LRT.inspectLg=false;
pars.LRT.subspacedim1=1;
pars.LRT.subspacedim2=1;
pars.LRT.G.precon=true;
pars.LRT.G.maxiter = 10;
pars.LRT.scaleksp=0;
pars.LRT.niter = 3;
LRT_nav_mask_enable = strcmp( pars.method, 'LRT');
lrt_nav_weight = 2e-2;
    
image_corrected = msDWIrecon(kspa, squeeze(TSE_sense_map), [], pars);  %no phase error maps needed here




%===================================================FINISH=================================================================













































%% APPENDEX Spiral NUFFT recon.
disp(' Spiral NUFFT recon...');
save_mat_fn = 'Sc08.mat';
close all;
[kx_length ch_nr shot_nr, dyn_nr] = size(nav_k_spa_data);

offcenter_xy = [0 0]; 
FOV_xy = [250 164.7727];
% nav_im_recon_nufft = [];
dyn_recon = 1:dyn_nr;
for d = 1:length(dyn_recon)
    tic
    dyn  = dyn_recon(d);
    disp(['dynamic: ',num2str(dyn)]);
    %=============== recon parameters =========================
    recon_par.ignore_kz = 0;
    recon_par.acq_dim = [42 42 13];  
    recon_par.recon_dim  = [42 42 13];
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




