%% load trajectry
clear; close all;
load('traj_Sc10_12_13.mat');
trj_meas_kx_t = trj_meas_kx(1:1222,1);
trj_meas_ky_t = trj_meas_ky(1:1222,1);
trj_meas_kz_t = trj_meas_kz(1:1222,1);

clear trj_meas_kx trj_meas_ky trj_meas_kz
trj_meas_kx = trj_meas_kx_t;
trj_meas_ky = trj_meas_ky_t;
trj_meas_kz = trj_meas_kz_t;




%% generate bart_phantom traj

xy_res = 30;
z_res = 15;
scale_foctor_kx_ky = max([(xy_res/2)./max(abs(trj_meas_kx)),(xy_res/2)./max(abs(trj_meas_ky))]);
scale_foctor_kz = (z_res/2)./max(abs(trj_meas_kz));
% scale_foctor_kx_ky = 18;
trj_meas_kx_scaled = trj_meas_kx * scale_foctor_kx_ky;
trj_meas_ky_scaled = trj_meas_ky * scale_foctor_kx_ky;
trj_meas_kz_scaled = trj_meas_kz * scale_foctor_kz;

traj_matrix_2D_phantom_rs(1,:) = trj_meas_kx_scaled;
traj_matrix_2D_phantom_rs(2,:) = trj_meas_ky_scaled;
traj_matrix_2D_phantom_rs(3,:) = zeros(1,length(trj_meas_kx_scaled));

traj_matrix_3D_phantom_rs(1,:) = trj_meas_kx_scaled;
traj_matrix_3D_phantom_rs(2,:) = trj_meas_ky_scaled;
traj_matrix_3D_phantom_rs(3,:) = trj_meas_kz_scaled;

%% generate bart_nufft traj

xy_res = 20;
z_res = 15;
scale_foctor_kx_ky = max([(xy_res/2)./max(abs(trj_meas_kx)),(xy_res/2)./max(abs(trj_meas_ky))]);
scale_foctor_kz = (z_res/2)./max(abs(trj_meas_kz));
% scale_foctor_kx_ky = 18;
trj_meas_kx_scaled = trj_meas_kx * scale_foctor_kx_ky;
trj_meas_ky_scaled = trj_meas_ky * scale_foctor_kx_ky;
trj_meas_kz_scaled = trj_meas_kz * scale_foctor_kz;

traj_matrix_2D_nufft_rs(1,:) = trj_meas_kx_scaled;
traj_matrix_2D_nufft_rs(2,:) = trj_meas_ky_scaled;
traj_matrix_2D_nufft_rs(3,:) = zeros(1,length(trj_meas_kx_scaled));

traj_matrix_3D_nufft_rs(1,:) = trj_meas_kx_scaled;
traj_matrix_3D_nufft_rs(2,:) = trj_meas_ky_scaled;
traj_matrix_3D_nufft_rs(3,:) = trj_meas_kz_scaled;

%% 2D simulation
xy_res = 30 
K_space_ideal = bart('phantom -s8 -k -t', traj_matrix_2D_phantom_rs);

bart_com_1 = sprintf('nufft -i -l0.01 -d%d:%d:1 ',xy_res,xy_res);
ima_nufft = bart(bart_com_1, traj_matrix_2D_nufft_rs, K_space_ideal);

figure(2001); montage(permute(abs(ima_nufft),[1 2 3 4]),'displayrange',[]); colormap gray

% estimate sense map
kspa_ideal_recon = bart('fft -u 7', ima_nufft);
sens_self_ecalib1 =bart('ecalib -m1', kspa_ideal_recon);

% get external sense map
k_cartesian = bart('phantom -s8 -x64 -k');

bart_com_2 = sprintf('resize -c 0 %d 1 %d',xy_res,xy_res);
k_cartesian_rs = bart(bart_com_2, k_cartesian);
sens_maps_cartesian = bart('ecalib -m1', k_cartesian_rs);

%sense unwrapping
ima_trj_pics = squeeze(bart('pics -S -l2 -r0.001 -t', ...
    traj_matrix_2D_nufft_rs, K_space_ideal, sens_maps_cartesian)); %or sens_maps_cartesian

figure(20021); 
subplot(121);imshow(abs(squeeze(bart('rss 8',ima_nufft))),[]); title('reference')
subplot(122);imshow(abs(ima_trj_pics),[]); title('pics+traj correction (spiral sense)')

%% 2D NUFFT simulation
%trajectory should scale to [-pi pi] = bart_nufft interval = 1; 
%scale to e.g. [-2pi 2pi] = bart_nufft interval = 2  >>> big FOV and image repetition


% cart_traj_nufft = pi/15*traj_matrix_2D_phantom_rs';
%or
cart_traj_nufft = pi/10*traj_matrix_2D_nufft_rs';


%correct direct recon
 A=nuFTOperator(cart_traj_nufft(:,1:2),[64, 64],ones(64, 64,8),6); 
 im_recon_nufft=regularizedReconstruction(A,col(K_space_ideal(1,:,:,:)),@L2Norm,0.5,'maxit',25);
       
 figure; imshow(abs(im_recon_nufft),[])      
 
 %correct sense recon
  A=nuFTOperator(cart_traj_nufft(:,1:2),[30, 30],squeeze(sens_maps_cartesian),6); 
 im_recon_nufft=regularizedReconstruction(A,col(K_space_ideal(1,:,:,:)),@L2Norm,0.5,'maxit',25);
       
 figure; imshow(abs(im_recon_nufft),[])    

%% 3D simulation

Kspa_spira3D = bart('phantom -s8 -k -3 -t', traj_matrix_3D_phantom_rs);

bart_com_3 = sprintf('nufft -i -l0.01 -d%d:%d:%d ',xy_res,xy_res,z_res);
ima_nufft_3D = bart(bart_com_3, traj_matrix_2D_phantom_rs, Kspa_spira3D);

figure(3001); montage(permute(abs(ima_nufft_3D(:,:,7,:)),[1 2 3 4]),'displayrange',[]); colormap gray

Kspa_spira3D_grid = bart('fft -u 7', ima_nufft_3D);
sens_3D_ecalib1 =bart('ecalib -c0 -m1', Kspa_spira3D_grid);


k_cartesian_3D = bart('phantom -s8 -x64 -3 -k');

bart_com_4 = sprintf('resize -c 0 %d 1 %d 2 %d',xy_res,xy_res,z_res);
k_cartesian_3D_rs = bart(bart_com_4, k_cartesian_3D);

sens_maps_3D_cartesian = bart('ecalib -c0 -m1', k_cartesian_3D_rs);

figure(3002); 
subplot(121)
montage(permute(abs(sens_3D_ecalib1(:,:,7,:)),[1 2 3 4]),'displayrange',[]); colormap gray
subplot(122);
montage(permute(abs(sens_maps_3D_cartesian(:,:,7,:)),[1 2 3 4]),'displayrange',[]); colormap gray

external_sense_recon = 1;

if(external_sense_recon)
    pic_recon_sense_map = sens_maps_3D_cartesian;
else
    pic_recon_sense_map = sens_3D_ecalib1;
end

ima_trj_pics = squeeze(bart('pics -S -l2 -r0.001 -t', ...
    traj_matrix_3D_phantom_rs, Kspa_spira3D, pic_recon_sense_map));

figure(30021);
subplot(131);
montage(permute(abs(bart('rss 8',ima_nufft_3D )),[1 2 4 3]),'displayrange',[]);
title('ref recon');
subplot(132);
montage(permute(abs(ima_trj_pics),[1 2 4 3]),'displayrange',[]);
title('self sens recon');
subplot(133);
montage(permute(abs(pic_recon_sense_map(:,:,:,1)),[1 2 4 3]),'displayrange',[]); colormap gray
title('self sens map channel 1');

%%