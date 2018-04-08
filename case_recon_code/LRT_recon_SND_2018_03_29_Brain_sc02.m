
clear; clc; close all
cd('/home/qzhang/lood_storage/divi/Users/qzhang/Data/2018_03_29_DISPENSE_LRT_Brain')


%% SET path for all the following steps
clear; close all; clc

data_fn = 'di_29032018_1706497_2_2_wip_sc2_3d_snd_brain_4bV4.raw';
sense_ref_fn = 'di_29032018_1706107_1000_5_wip_senserefscanV4.raw';
coil_survey_fn  = 'di_29032018_1702382_1000_2_wip_coilsurveyscanV4.raw';





%% TSE imaging data loading (default recon)
close all; clc
disp(' TSE data sorting and default recon...')

parameter2read.dyn = [];

[ima_k_spa_data,TSE.ky_matched,TSE.kz_matched,TSE.shot_matched, TSE.ch_dim,ima_kspa_sorted, ima_default_recon, TSE_sense_map, TSE.kxrange, TSE.kyrange, TSE.kzrange, TSE.VirtualCoilMartix] = ...
    TSE_data_sortting(data_fn, sense_ref_fn, coil_survey_fn,parameter2read);    
    %use TSE_data_sortting_no_Nav when no navigator acquired

figure(610); immontage4D(permute(abs(ima_default_recon(:,:,:,:)),[1 2 4 3]), [0 200]);

TSE
assert(length(TSE.ky_matched)==size(ima_k_spa_data,2),'Profile number does not match with data size!')
disp('-finished- ');


%% SET parameter

%-------------TSE pars.-------------%
TSE.kxrange = [-352 -1];
TSE.dyn_dim = size(ima_kspa_sorted, 5);
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
pars.recon_x_locs = 100:380; %80:270;
pars.enabled_ch = 1:TSE.ch_dim;
pars.b0_shots = 1:50; % means first dynamic
pars.nonb0_shots = 51:100;
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
figure(3); montage(permute(squeeze(abs(TSE_sense_map(250,:,:,:))),[1 2 4 3]),'displayrange',[]); xlabel('SENSE maps')
%-------------------end---------------%

%% Get b0 data

%=======TSE imaging=======%
x_loc = 150;
hybrid_k_spa_data = zeros(TSE.kx_dim, size(ima_k_spa_data,2));
pad_left = floor((TSE.kx_dim - size(ima_k_spa_data,1))/2);
hybrid_k_spa_data(pad_left+1:pad_left+size(ima_k_spa_data,1),:) = ifft1d(ima_k_spa_data);

b0_kpa = sort_k_spa_sh_by_sh_3(squeeze(hybrid_k_spa_data(x_loc,:)), pars.b0_shots, TSE, pars);
%----combine shots directly; non-zero average in the shot dim
b0_kspa_comb = sum(b0_kpa, 5)./ sum(abs(b0_kpa)>0, 5); b0_kspa_comb(find(isnan(b0_kspa_comb)+isinf(b0_kspa_comb))) = 0;
b0_ima = bart('pics -RT:7:0:0.001 -d5', b0_kspa_comb, TSE_sense_map(x_loc,:,:,:));
figure(2); imshow(abs(squeeze(b0_ima)),[]);% xlabel('CS recon of b0 data: shot directly combined')
b0_kspa_full = squeeze(fft2d(bsxfun(@times, b0_ima, TSE_sense_map)));

size(b0_kspa_full) 
% clear b0_kpa b0_kspa b0_ima



%% Get non-b0 data

%=======TSE imaging=======%
nonb0_kpa = sort_k_spa_sh_by_sh_3(squeeze(hybrid_k_spa_data(x_loc,:)), pars.nonb0_shots, TSE, pars);


size(nonb0_kpa) 






%% RECON

% kspa = squeeze(cat(5, b0_kspa, nonb0_kpa));
kspa = squeeze(cat(5, b0_kspa));

kspa = repmat(kspa, [1 1 1 1 2]);

%------LRT pars---------%
pars.method='LRT'; %POCS_ICE CG_SENSE_I CG_SENSE_K LRT

pars.LRT.Lg=3;
pars.LRT.L3=3;
pars.LRT.L4=1;
pars.LRT.mu = 1e3;
pars.LRT.beta = 1;
pars.LRT.lambda = 9e-3;

pars.LRT.sparsity_transform='TV';
%     pars.LRT.Imref=cat(3, repmat(squeeze(im_b0_ref(recon_x_loc,:,:,:)), [1 1 1 2]), repmat(squeeze(im_ref(recon_x_loc,:,:,:)), [1 1 length(recon_shot_range) 2]));
pars.LRT.x=30;  %loc in dim2
pars.LRT.y=90; %loc in dim1
pars.LRT.increase_penalty_parameters=false;
pars.LRT.inspectLg=false;
pars.LRT.subspacedim1=1;
pars.LRT.subspacedim2=1;
pars.LRT.G.precon=true;
pars.LRT.G.maxiter = 10;
pars.LRT.scaleksp=0;
pars.LRT.niter = 5;

    
image_corrected = msDWIrecon(kspa, squeeze(TSE_sense_map(x_loc,:,:,:)), [], pars);  %no phase error maps needed here

figure(673); subplot(221); imshow((squeeze(abs(image_corrected(:,:,:,1,1)))),[]); title('b0 nav')
subplot(222); imshow((squeeze(abs(image_corrected(:,:,:,1,2)))),[]); title('b0 TSE')
figure(673); subplot(223); imshow((squeeze(angle(image_corrected(:,:,:,1,1)))),[]); title('phase b0 nav')
subplot(224); imshow((squeeze(angle(image_corrected(:,:,:,1,2)))),[]); title('phase b0 TSE')

figure(674); montage(permute(squeeze(abs(image_corrected(:,:,:,2:end,1))),[1 2 4 3]),'displayrange',[]); title('nav')
figure(675); montage(permute(squeeze(abs(image_corrected(:,:,:,2:end,2))),[1 2 4 3]),'displayrange',[]); title('TSE')

figure(676); montage(permute(squeeze(angle(image_corrected(:,:,:,2:end,1))),[1 2 4 3]),'displayrange',[]); colormap jet; title('phase nav')
figure(677); montage(permute(squeeze(angle(image_corrected(:,:,:,2:end,2))),[1 2 4 3]),'displayrange',[]); colormap jet; title('phase TSE')

% pars.LRT.NUFFT_nav_sense = pars.LRT.NUFFT_nav_1ch;
% image_corrected = msDWIrecon(kspa_1ch, squeeze(ones(size(TSE_sense_map(:,:,1,1)))), [], pars);  %no phase error maps needed here


%==========================================FINISH===============================================================
