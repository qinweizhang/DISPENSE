%gen 2D traj
id_2d = find(traj_matrix_ideal_rs(3,:,:) == 0);
temp = traj_matrix_ideal_rs(1,:,:);
traj_matrix_ideal_rs_2D(1,:) =temp(id_2d);
temp = traj_matrix_ideal_rs(2,:,:);
traj_matrix_ideal_rs_2D(2,:) =temp(id_2d);
traj_matrix_ideal_rs_2D(3,:) =zeros(1,length(id_2d));

id_2d = find(traj_matrix_cor_rs(3,:,:) == 0);
temp = traj_matrix_cor_rs(1,:,:);
traj_matrix_shifted_rs_2D(1,:) =temp(id_2d);
temp = traj_matrix_cor_rs(2,:,:);
traj_matrix_shifted_rs_2D(2,:) =temp(id_2d);
traj_matrix_shifted_rs_2D(3,:) =zeros(1,length(id_2d));

%% 2D simulation
% ideal
K_space_ideal = bart('phantom -k -s 8 -t', traj_matrix_ideal_rs_2D);

ima_ideal = bart('nufft -i -l0.01 ', traj_matrix_ideal_rs_2D, K_space_ideal);
% ima_ideal_fftshifted = fftshift(fftshift(ima_ideal,2),3);
kspa_ideal_recon = bart('fft -u 7', ima_ideal);
figure(2001); montage(permute(abs(ima_ideal),[1 2 3 4]),'displayrange',[]); colormap gray
sens_ecalib =bart('ecalib -m1', kspa_ideal_recon);

% shifted
traj_matrix_shifted_rs_2D_tailored = traj_matrix_shifted_rs_2D(:,1:2:end);
K_space_trj_shifted = bart('phantom -k -s 8 -t', traj_matrix_shifted_rs_2D_tailored);

ima_trj_shifted = bart('nufft -i -l0.01 ', traj_matrix_shifted_rs_2D_tailored, K_space_trj_shifted);
kspa_trj_shifted_recon = bart('fft -u 7', ima_trj_shifted);

figure(2002); montage(permute(abs(ima_trj_shifted),[1 2 3 4]),'displayrange',[]); colormap gray


% sens_ecalib_rs = permute(sens_ecalib, [5 1 2 3 4]);
ima_trj_shifted_pics = squeeze(bart('pics -S -l2  -r0.01 -t', traj_matrix_shifted_rs_2D_tailored, K_space_trj_shifted, sens_ecalib));
figure(20021); 
subplot(131);imshow(abs(squeeze(bart('rss 8',ima_ideal))),[]); title('reference')
subplot(132);imshow(abs(squeeze(bart('rss 8',ima_trj_shifted))),[]); title('uncorrection')
subplot(133);imshow(abs(ima_trj_shifted_pics),[]); title('pics+traj correction')

%% 3D simulation
% ideal
K_space_ideal = bart('phantom -k -3 -s 8 -t', traj_matrix_ideal_rs);

ima_ideal = bart('nufft -i -l0.01 ', traj_matrix_ideal_rs, K_space_ideal);
% ima_ideal_fftshifted = fftshift(fftshift(ima_ideal,2),3);

figure(3001); montage(permute(abs(ima_ideal(:,:,1,:)),[1 2 3 4]),'displayrange',[]); colormap gray
for z = 1: size(kspa_ideal_recon,3)
    kspa_ideal_recon = bart('fft -u 7', ima_ideal(:,:,z,:));
    sens_ecalib(:,:,z,:) =bart('ecalib -m1', kspa_ideal_recon);
end
figure(30011); montage(permute(abs(sens_ecalib(:,:,2,:)),[1 2 3 4]),'displayrange',[]); colormap gray


% shifted
traj_matrix_cor_rs_tailored = traj_matrix_cor_rs(:,1:1:end);
K_space_trj_shifted = bart('phantom -k -3 -s 8 -t', traj_matrix_cor_rs_tailored);

ima_trj_shifted = bart('nufft -i -l0.1 ', traj_matrix_cor_rs_tailored, K_space_trj_shifted);
kspa_trj_shifted_recon = bart('fft -u 7', ima_trj_shifted);

figure(3002); montage(permute(abs(ima_trj_shifted(:,:,2,:)),[1 2 3 4]),'displayrange',[]); colormap gray


% sens_ecalib_rs = permute(sens_ecalib, [5 1 2 3 4]);
ima_trj_shifted_pics = squeeze(bart('pics -S -l2  -r1 -t', traj_matrix_cor_rs_tailored, K_space_trj_shifted, sens_ecalib));
figure(30021); 
subplot(131);imagesc(array2mosaic(abs(squeeze(bart('rss 8',ima_ideal))))); title('reference'); colormap gray; axis equal; axis off
subplot(132);imagesc(array2mosaic(abs(squeeze(bart('rss 8',ima_trj_shifted))))); title('reference'); colormap gray; axis equal; axis off
subplot(133);imagesc(array2mosaic(abs(ima_trj_shifted_pics))); title('reference'); colormap gray; axis equal; axis off


%% display
disp_kz = 0;
figure(2003);

disp_range = find(traj_matrix_ideal_rs(3,:,:) == disp_kz);

temp_x = traj_matrix_ideal_rs(1,disp_range); temp_y = traj_matrix_ideal_rs(2,disp_range); temp_k = abs(K_space_ideal(disp_range));
subplot(121); scatter(temp_x(:), temp_y(:),30, temp_k(:),'filled'); title('uncor mag'); xlabel('kx'); ylabel('ky')
xlim([-80 80]); ylim([-80 80]);
disp_range = find(traj_matrix_cor_rs(3,:,:) == disp_kz);
temp_x_shifted = traj_matrix_cor_rs(1,disp_range); temp_y_shifted = traj_matrix_cor_rs(2,disp_range); temp_k_shifted = abs(K_space_trj_shifted(disp_range));
subplot(122); scatter(temp_x_shifted(:), temp_y_shifted(:),30, temp_k_shifted(:),'filled');title('cor mag'); xlabel('kx'); ylabel('ky')
xlim([-80 80]); ylim([-80 80]);

figure(2004);

disp_range = find(traj_matrix_ideal_rs(3,:,:) == disp_kz);

temp_x = traj_matrix_ideal_rs(1,disp_range); temp_y = traj_matrix_ideal_rs(2,disp_range); temp_k = angle(K_space_ideal(disp_range));
subplot(121); scatter(temp_x(:), temp_y(:),30, temp_k(:),'filled'); title('uncor mag'); xlabel('kx'); ylabel('ky')
xlim([-80 80]); ylim([-80 80]);
disp_range = find(traj_matrix_cor_rs(3,:,:) == disp_kz);
temp_x_shifted = traj_matrix_cor_rs(1,disp_range); temp_y_shifted = traj_matrix_cor_rs(2,disp_range); temp_k_shifted = angle(K_space_trj_shifted(disp_range));
subplot(122); scatter(temp_x_shifted(:), temp_y_shifted(:),30, temp_k_shifted(:),'filled');title('cor mag'); xlabel('kx'); ylabel('ky')
xlim([-80 80]); ylim([-80 80]);

figure(2005);
subplot(121); imshow(squeeze(abs(kspa_trj_shifted_recon(:,:,2)))',[]);title('nufft recon mag'); xlabel('kx'); ylabel('ky'); colormap jet
subplot(122); imshow(squeeze(angle(kspa_ideal_recon(:,:,2)))',[]);title('nufft recon mag'); xlabel('kx'); ylabel('ky'); colormap jet

%%
traj_nufft_cor_actual = traj_nufft_cor;
% ideal
K_space_data = bart('phantom -k -3 -s 8 -t', traj_nufft_cor_actual);
ima_trj_shifted = bart('nufft -i -l0.1 ', traj_nufft_cor_actual, K_space_data);
ima_trj_shifted_rss = bart('rss 8',ima_trj_shifted);
figure(10); montage(permute(abs(ima_trj_shifted(:,:,2,:)),[1 2 3 4]),'displayrange',[])
figure(11); montage(permute(abs(ima_trj_shifted_rss),[1 2 4 3]),'displayrange',[])

% nufft
% scale trajectory [-pi pi]
nufft_dim = round(max(traj_nufft_cor,[],2) - min(traj_nufft_cor,[],2)+1); 
nufft_dim(3) = 2;

scale_factor_x = pi / max(abs(traj_nufft_cor(1,:)));
scale_factor_y = pi / max(abs(traj_nufft_cor(2,:)));
scale_factor_z = pi / max(abs(traj_nufft_cor(3,:)));

traj_nufft_cor_scaled(:,1) = traj_nufft_cor(1,:) * scale_factor_x;
traj_nufft_cor_scaled(:,2) = traj_nufft_cor(2,:) * scale_factor_y;
traj_nufft_cor_scaled(:,3) = traj_nufft_cor(3,:) * scale_factor_z;

clear im_recon_nufft_cor im_recon_nufft_cor_fftshifted


%     A=nuFTOperator(traj_nufft_cor_scaled,[nufft_dim(1), nufft_dim(2), nufft_dim(3)],ones(nufft_dim(1), nufft_dim(2), nufft_dim(3),ch_dim),oversample_factor);
A=nuFTOperator(traj_nufft_cor_scaled(:,:),[nufft_dim(1), nufft_dim(2), nufft_dim(3)],ones(nufft_dim(1), nufft_dim(2), nufft_dim(3)),oversample_factor);




im_recon_nufft_cor=regularizedReconstruction(A,col(K_space_data),@L2Norm,0.5,'maxit',25);
% im_recon_nufft_cor=A'*col(K_space_data);

im_recon_nufft_cor_fftshifted = fftshift(fftshift(im_recon_nufft_cor,2),3);
figure(25); montage(permute(abs(im_recon_nufft_cor ),[1 2 4 3]),'displayrange',[]); title('cor')



