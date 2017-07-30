clear; close all; clc;
fn = 'data_Sc19_2D_cor.mat';
 load('traj_Sc8-9.mat'); load(fn);

%scle trajectory to match kx_range, ky_range (-pi to pi)

shot_nr = 1; 
diffusion_nr = 1;
nav_k_spa_data = squeeze(nav_k_spa_data);
skip_point = 0;
end_point = 3000%length(trj_meas_kx);
selected_point = [skip_point+1:end_point];

trj_meas_kx_t = squeeze(trj_meas_kx(selected_point,1));
trj_meas_ky_t = squeeze(trj_meas_ky(selected_point,1));
trj_meas_kz_t = squeeze(trj_meas_kz(selected_point,1));


clear im_recon_nufft
for shot_nr = 1: 500%size(nav_k_spa_data,3)
% for shot_nr = 1: 3
%     shot_nr = 1 + (shot_nr - 1) *5; 
    sig_kspa = nav_k_spa_data(selected_point,:,shot_nr);
    
    %trajectory correction
    clear trj_meas_kx trj_meas_ky trj_meas_kz
    trj_meas_kx = trj_meas_kx_t;
    trj_meas_ky = trj_meas_ky_t;
    trj_meas_kz = trj_meas_kz_t;

    scale_foctor = max(pi/max(abs(trj_meas_kx(:,diffusion_nr))),pi/max(abs(trj_meas_ky(:,diffusion_nr))));
    trj_meas_kx_scaled = trj_meas_kx(:,diffusion_nr) * scale_foctor;
    trj_meas_ky_scaled = trj_meas_ky(:,diffusion_nr) * scale_foctor;

    trj_meas_kx_scaled = trj_meas_kx_scaled - 2 * pi / 24 * linear_phase_error_input_x_in_TSEksp_pixel(shot_nr );
    trj_meas_ky_scaled = trj_meas_ky_scaled - 2 * pi / 24 * linear_phase_error_input_y_in_TSEksp_pixel(shot_nr );
    

    trj_nufft = double(cat(2, trj_meas_kx_scaled, trj_meas_ky_scaled));


    disp(['shot nr: ',num2str(shot_nr)]);
    for ch =1:size(nav_k_spa_data,2)
        sig_nufft = double(sig_kspa(:, ch));
        sig_nufft = sig_nufft';
        A=nuFTOperator(trj_nufft,[24, 24],ones(24, 24),6);

        % simple inverse
%          temp = A'*sig_nufft';



         %call CG-SENSE with L2-norm regularization
         temp = regularizedReconstruction(A,sig_nufft',@L2Norm,0.5,'maxit',25);
         nav_im_recon_nufft(:,:,ch,shot_nr) = temp; 
    %     im_recon_nufft(:,:,ch)=regularizedReconstruction(A,sig_nufft','maxit',25);
    end
    figure(30);  imagesc(angle(temp), [-pi pi]); axis off; axis equal; colormap jet;
end
ch_nr = 6;
figure(35); montage(abs(nav_im_recon_nufft(:,:,ch_nr,:)),'displayrange',[]); title(['all shots from ch', num2str(ch_nr)]);
figure(36); montage(angle(nav_im_recon_nufft(:,:,ch_nr,:)),'displayrange',[-pi pi]); colormap jet; title(['all shots from ch', num2str(ch_nr)]);

phase_test_ch = (angle(nav_im_recon_nufft(:,:,ch_nr,:)));
for shot_nr = 1:size(nav_k_spa_data,3)
    phase_test_ch_diff(:,:,1,shot_nr) = phase_test_ch(:,:,1,shot_nr) -phase_test_ch(:,:,1,4) ;
end
figure(37); montage(phase_test_ch_diff,'displayrange',[-2*pi 2*pi]); colormap jet; title(['all shots phase differece from ch', num2str(ch_nr)]);

nav_im_recon_nufft_rss = bart('rss 4',nav_im_recon_nufft);
figure(31);montage(abs(nav_im_recon_nufft_rss),'Displayrange',[]); title(['all shots rss']);
% figure(32);montage(angle(permute(im_recon_nufft,[1 2 4 3])),'Displayrange',[-pi pi]); colormap 'jet'
% figure(33); imshow(abs(((im_recon_nufft_rss))),[])

save(fn,'nav_im_recon_nufft','nav_im_recon_nufft_rss','-append');