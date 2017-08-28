%     Use msDWIrecon function to correct phase inherence during multishot DPsti-TSE
%
%     INPUT
%
%     ima_k_spa_data:             Unsorted TSE kspace data in [kx, nprofiles]
%     TSE:                        Structure contains labels (ky, kz, shot, channel, dyn) for every profiles in ima_k_spa_data
%     TSE_sense_map:              sense maps from sense reference scan
%     nav_im:                     navigator cpx images for every shot in [nav_x nav_y nav_z shots]
%
%     OUTPUT
%
%     image_corrected:            Corrected DPsti-TSE images
%
%     (c) Qinwei Zhang (q.zhang@amc.uva.nl) 2017 @AMC Amsterdam

function image_corrected = DPsti_TSE_phase_error_cor(ima_k_spa_data, TSE, TSE_sense_map, nav_im, pars)
%% Data check
max_shot = max(TSE.shot_matched);
assert(max_shot == size(nav_im, 4));
assert(TSE.ch_dim == size(TSE_sense_map, 4)||isempty(TSE_sense_map));

profiles_per_dyn = size(ima_k_spa_data,2)/TSE.dyn_dim;
assert(round(profiles_per_dyn) == profiles_per_dyn);
shots_per_dyn = max_shot/TSE.dyn_dim;
assert(round(shots_per_dyn) == shots_per_dyn);
%% Preprocessing on kspace data for b=0
disp('Preprocessing on kspace data for b=0');
if(isempty(pars.b0_shots))
    b0_shots_range = 1:shots_per_dyn; %by default, the first dynamic
else
    b0_shots_range = pars.b0_shots;
end

tic;
kx_dim = TSE.kxrange(2) - TSE.kxrange(1) + 1;
assert(kx_dim ==  size(ima_k_spa_data,1));
max_ky = max(abs(TSE.kyrange)); 
max_kz = max(abs(TSE.kzrange));
ky_dim = TSE.kyrange(2) - TSE.kyrange(1) + 1;  %max_ky * 2 + 1; 
kz_dim = TSE.kzrange(2) - TSE.kzrange(1) + 1;
ch_dim = TSE.ch_dim;
sh_dim = range(TSE.shot_matched)+1;
kspa_all = zeros(kx_dim, ky_dim, kz_dim, ch_dim, sh_dim);

for prof_idx = 1:size(ima_k_spa_data, 2)
    ky_idx = TSE.ky_matched(prof_idx) - TSE.kyrange(1) + 1;
    kz_idx = TSE.kz_matched(prof_idx) - TSE.kzrange(1) + 1;
    ch_idx = mod(prof_idx, ch_dim) + 1;
    sh_idx = TSE.shot_matched(prof_idx);
    kspa_all(:,ky_idx,kz_idx,ch_idx, sh_idx) = ima_k_spa_data(:,prof_idx);
end
size(kspa_all)

kspa_all = kspa_all(:,:,:,pars.enabled_ch,:);

% remove stupid checkerboard pattern
che=create_checkerboard([1,size(kspa_all,2),size(kspa_all,3)]);
kspa_all=bsxfun(@times,kspa_all,che);
kspa_all=squeeze(kspa_all);

kspa_b0 = sum(kspa_all(:,:,:,:,b0_shots_range), 5); %4D b0 kspace [kx ky kz nc]

pp=bart('fft -i 7',kspa_b0);
im_b0=bart('rss 8',pp);
figure(1); montage(permute(abs(im_b0),[1 2 4 3]),'displayrange',[]); title('b0 images');
clear pp
toc;

%% Preprocessing on sense data
disp('Preprocessing on sense data');
tic

%estimate sense maps

if(strcmp(pars.sense_map, 'ecalib'))
    
    ecalib_sense_map_3D = bart('ecalib -S -m1', kspa_b0);
    figure(3);  montage(angle(ecalib_sense_map_3D(:,:,8,:)),'displayrange',[]); title('ecalib sense map')
    figure(4);  montage(abs(ecalib_sense_map_3D(:,:,8,:)),'displayrange',[]); title('ecalib sense map')
    ecalib_sense_map_3D = normalize_sense_map(ecalib_sense_map_3D);
    
    sense_map_3D = ecalib_sense_map_3D;
    
    
    %match the size of TSE_sens_map to kspa_xyz
    %TODO
elseif(strcmp(pars.sense_map, 'external'))
    if(isempty(TSE_sense_map))
        dim = [size(kspa_b0, 2) size(kspa_b0, 2) size(kspa_b0,3)];
        os = [1, 1, 1];
        TSE_sense_map = get_sense_map_external(pars.sense_ref, pars.data_fn, pars.coil_survey, dim, os);
        rs_command = sprintf('resize -c 0 %d', size(kspa_b0, 1));
        TSE_sense_map = bart(rs_command, TSE_sense_map);
    end
    sense_map_3D = normalize_sense_map(TSE_sense_map(:,:,:,pars.enabled_ch ));
else
    error('sense map source not indentified.')
end
toc

%% Preprocessing on kspace data for nonb0
disp('Preprocessing on kspace data for nonb0');

% correction_shots_range = [shots_per_dyn+1: 2*shots_per_dyn];
if(isempty(pars.nonb0_shots))
    nonb0_shots_range = 1: shots_per_dyn;
else
    nonb0_shots_range = pars.nonb0_shots;
end


tic
kspa_xyz = kspa_all(:,:,:,:,nonb0_shots_range);

kk = sum(kspa_xyz, 5); %4D b0 kspace [kx ky kz nc]
pp=bart('fft -i 7',kk);
im_nonb0=bart('rss 8',pp);
figure(2); montage(permute(abs(im_nonb0),[1 2 4 3]),'displayrange',[]); title('direct recon');
clear kk pp

%to hybrid space
clear kspa_x_yz;
kspa_x_yz = ifft1d(kspa_xyz);
toc



%% preprocssing on phase error data
disp('Preprocessing on phase error data');
tic
%match size of nav_im and to kspa_xyz
nav_im_1 = nav_im(:,:,:,nonb0_shots_range);
[kx, ky, kz, nc, nshot] = size(kspa_xyz);
assert(size(nav_im_1, 4)==nshot);

nav_k_1 = bart('fft 7', nav_im_1);
resize_command = sprintf('resize -c 0 %d 1 %d 2 %d 3 %d', ky, ky, kz, nshot); %nav and TSE have the same FOV, but TSE have oversampling in x, so use ky instead of kx for the 1st dimension
nav_k_1 = bart(resize_command, nav_k_1);
nav_im_2 = bart('fft -i 7', nav_k_1);

resize_command_2 = sprintf('resize -c 0 %d 1 %d 2 %d 3 %d', kx, ky, kz, nshot);
nav_im_2 = bart(resize_command_2, nav_im_2);
phase_error_3D = nav_im_2;
phase_error_3D = normalize_sense_map(phase_error_3D); %miss use normalize_sense_map

figure(5);
immontage4D(angle(phase_error_3D),[-pi pi]); colormap jet; title('phase error maps int.')
xlabel('shot #'); ylabel('slice locations');


toc
%% msDWI recon
disp('recon');
image_corrected = zeros(kx, ky, kz);
tic;
recon_x = [1:128]; %pars.recon_x_locs;
for x_idx = 1:length(recon_x)
    
    recon_x_loc = recon_x(x_idx);
    %recon parameter
    pars.lamda=0.1;
    pars.nit=10;
    pars.method='CG_SENSE';
    
    %=========select data. fixed=============================================
    kspa = squeeze(kspa_x_yz(recon_x_loc, :, :, :, :));
    sense_map = squeeze(sense_map_3D(recon_x_loc,:,:,:));
    phase_error = permute(squeeze(phase_error_3D(recon_x_loc,:,:,:,:)),[1 2 4 3]);
    %========================================================================
    
    image_corrected(recon_x_loc,:,:) = msDWIrecon(kspa, sense_map, phase_error, pars);
    
    image_corrected(isnan(image_corrected)) = 0;
    
    %display
    figure(101);
    subplot(131);imshow(squeeze(abs(im_b0(recon_x_loc,:,:))),[]); title('b0');
    subplot(132);imshow(squeeze(abs(im_nonb0(recon_x_loc,:,:))),[]); title('direct recon');
    subplot(133);imshow(squeeze(abs(image_corrected(recon_x_loc,:,:))),[]); title('msDWIrecon');
    
end

toc;

    %display
    figure(102);
    subplot(131); montage(permute(abs(im_b0),[1 2 4 3]),'displayrange',[]); title('b0');
    subplot(132); montage(permute(abs(im_nonb0),[1 2 4 3]),'displayrange',[]); title('direct recon');
    subplot(133); montage(permute(abs(image_corrected),[1 2 4 3]),'displayrange',[]); title('msDWIrecon');
end