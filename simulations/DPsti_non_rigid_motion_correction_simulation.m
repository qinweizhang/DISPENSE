%==========DPsti-TSE nonrigid motion correction simulation============
cd('/home/qzhang/lood_storage/divi/Users/qzhang/SoSND/SoSND/simulations')
clear; close all; clc
load('test_data.mat')

%% reference
kx_dim = size(ima_k_spa_ideal,1);
ky_dim = max(TSE.ky_matched) * 2 + 1; %consider half scan
kz_dim = max(TSE.kz_matched) * 2 + 1; %consider half scan
ch_dim = TSE.ch_dim; 
kspa_xyz_ref = zeros(kx_dim, ky_dim, kz_dim, TSE.ch_dim);

for prof_idx = 1:size(ima_k_spa_ideal, 2)
    ky_idx = TSE.ky_matched(prof_idx) + max(TSE.ky_matched) + 1;
    kz_idx = TSE.kz_matched(prof_idx) + max(TSE.kz_matched) + 1;
    ch_idx = mod(prof_idx, ch_dim) + 1;
    kspa_xyz_ref(:,ky_idx,kz_idx,ch_idx) = ima_k_spa_ideal(:,prof_idx); 
end
size(kspa_xyz_ref)

% remove stupid checkerboard pattern
che=create_checkerboard([1,size(kspa_xyz_ref,2),size(kspa_xyz_ref,3)]);
kspa_xyz_ref=bsxfun(@times,kspa_xyz_ref,che);
kspa_xyz_ref=squeeze(kspa_xyz_ref);

ima_ref = bart('fft -i 7',kspa_xyz_ref);
ima_coil_combined = bart('rss 8', ima_ref);
figure(1); montage(permute(abs(ima_coil_combined),[1 2 4 3]),'displayrange',[]);


sense_map = bart('ecalib -S -m1', kspa_xyz_ref);
figure(2);  montage(angle(sense_map(:,:,25,:)),'displayrange',[]);
figure(3);  montage(abs(sense_map(:,:,25,:)),'displayrange',[]);
sense_map = normalize_sense_map(sense_map);

clear TSE_sens_map ima_coil_combined

%% simulated phase error data
kx_dim = size(ima_k_spa_ideal,1);
ky_dim = max(TSE.ky_matched) * 2 + 1; %consider half scan
kz_dim = max(TSE.kz_matched) * 2 + 1; %consider half scan
ch_dim = TSE.ch_dim; 
sh_dim = range(TSE.shot_matched)+1;
kspa_xyz = zeros(kx_dim, ky_dim, kz_dim, ch_dim, sh_dim);

for prof_idx = 1:size(ima_k_spa_ideal, 2)
    ky_idx = TSE.ky_matched(prof_idx) + max(TSE.ky_matched) + 1;
    kz_idx = TSE.kz_matched(prof_idx) + max(TSE.kz_matched) + 1;
    ch_idx = mod(prof_idx, ch_dim) + 1;
    sh_idx = TSE.shot_matched(prof_idx);
    kspa_xyz(:,ky_idx,kz_idx,ch_idx, sh_idx) = ima_k_spa_ideal(:,prof_idx); 
end
size(kspa_xyz)

% remove stupid checkerboard pattern
che=create_checkerboard([1,size(kspa_xyz,2),size(kspa_xyz,3)]);
kspa_xyz=bsxfun(@times,kspa_xyz,che);
kspa_xyz=squeeze(kspa_xyz);

%phase error for every shot
phase_error = exp(i .* zeros(kx_dim, ky_dim, kz_dim, 1, sh_dim));

%corrupted kspa_xyz
kspa_xyz = bsxfun(@times,kspa_xyz,phase_error);

%direct recon
kk=sum(kspa_xyz,5);
pp=bart('fft -i 7',kk);
im_recon_direct=bart('rss 8',pp);
figure(2); montage(permute(abs(im_recon_direct),[1 2 4 3]),'displayrange',[]); title('direct recon');
clear kk pp im_recon_direct
%% recon
recon_x_loc = 100;
kspa_xyz = ifft1d(kspa_xyz);
kspa = squeeze(kspa_xyz(recon_x_loc, :, :, :, :));
sense_map = squeeze(sense_map_3D(recon_x_loc,:,:,:));
phase_error = permute(squeeze(phase_error_3D(recon_x_loc,:,:,:,:)),[1 2 4 3]);

image_corrected = msDWIrecon(kspa, sense_map, phase_error);