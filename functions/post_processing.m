%% set Matlab path and TOOLBOX_PATH environment variable
clear; close all;
clc
current_mat_file = 'data_Sc10_2D.mat';
raw_data_fn = 'dp_17052017_1418250_10_2_wipsc23ddpnavlinearexperiment1senseV4.raw';
%% 3D sprial navigator data phase unwrapping
% optional not part of the recon pipeline

% =================== load image
load(current_mat_file)
[nav_dimx, nav_dimy, nav_chan, nav_shot] = size(nav_im_recon_nufft);



%------------------------------------------------------PHSE UNWRAPPING----------------------------------%

Manually_seed_point = 0;
for ch_idx = 1:nav_chan
    for shot_idx = 1:nav_shot
        nav_ima_phase_unwrapped(:,:,ch_idx,shot_idx) = GoldsteinUnwrap2D_r1_func(nav_im_recon_nufft(:,:,ch_idx,shot_idx), Manually_seed_point);
        figure(22); 
        subplot(121); imshow(squeeze(angle(nav_im_recon_nufft(:,:,ch_idx,shot_idx))),[-pi pi]); title('L: wraped R: unwrapped ') ;  colorbar
        subplot(122); imshow(squeeze(nav_ima_phase_unwrapped(:,:,ch_idx,shot_idx)),[-pi pi]); colormap jet;title(['shot: ',num2str(shot_idx),'ch: ',num2str(ch_idx)]);  colorbar
    end
end



ch_idx = 6;
% ===================display wrapped phase
phase_wrapped = angle(nav_im_recon_nufft);

figure(601);

montage((phase_wrapped(:,:,ch_idx,:)),'displayrange',[-pi pi]); colormap jet; title(['all shots from ch', num2str(ch_idx)]);
% ===================display wrapped phase difference
phase_test_ch = phase_wrapped(:,:,ch_idx,:);
for shot_idx = 1:nav_shot
    phase_test_ch_diff(:,:,1,shot_idx) = phase_test_ch(:,:,1,shot_idx) -phase_test_ch(:,:,1,4) ;
end
figure(602); montage(phase_test_ch_diff,'displayrange',[-pi pi]); colormap jet; title(['all shots phase differece from ch', num2str(ch_idx)]);


% ====================display unwrapped phase difference
phase_test_ch = nav_ima_phase_unwrapped(:,:,ch_idx,:);
for shot_idx = 1:nav_shot
    phase_test_ch_diff(:,:,1,shot_idx) = phase_test_ch(:,:,1,shot_idx) -phase_test_ch(:,:,1,4) ;
end
figure(604); montage(phase_test_ch_diff,'displayrange',[-5*pi 5*pi]); colormap jet; title(['unwrapped phase difference:all shots phase differece from ch', num2str(ch_idx)]);

save(current_mat_file,'nav_ima_phase_unwrapped','-append');


%% DPsti-Data


MR_DPstiTSE = MRecon(raw_data_fn);

MR_DPstiTSE = MR_DPstiTSE.Copy;

MR_DPstiTSE.Parameter.Parameter2Read.typ = 1;
MR_DPstiTSE.Parameter.Parameter2Read.mix = 0;  
MR_DPstiTSE.Parameter.Recon.ImmediateAveraging = 'No';
MR_DPstiTSE.ReadData;
MR_DPstiTSE.RandomPhaseCorrection;
% MR_DPstiTSE.RemoveOversampling;
MR_DPstiTSE.PDACorrection;
MR_DPstiTSE.DcOffsetCorrection;
MR_DPstiTSE.MeasPhaseCorrection;

%--------------creat a label for shot number--------------&
    dyn=MR_DPstiTSE.Parameter.Labels.Index.dyn;
    typ=MR_DPstiTSE.Parameter.Labels.Index.typ;
    mix=MR_DPstiTSE.Parameter.Labels.Index.mix; 
    readdata_ix=find((typ==1)&(mix == 0));
    type_1_ix = find(typ==1);
    shot = zeros(length(typ),1);

    rf=MR_DPstiTSE.Parameter.Labels.Index.rf;
    first_echoes_ix = find(rf==0);
    ch_num = MR_DPstiTSE.Parameter.Labels.CoilNrs; %total number of the channels
    first_echoes_ty1_ix = first_echoes_ix(find(first_echoes_ix>=type_1_ix(1))); %first echoes for every shots
    first_echoes_ty1_first_ch = first_echoes_ty1_ix(1:length(ch_num)*2:end); %first lables for every shots; *2 means nav data 14ch + ima data 14ch
%     first_echoes_ty1_first_ch = first_echoes_ty1_ix(1:length(ch_num)*1:end); %first lables for every shots; when no DPnav enabled

    k_line_per_shot = first_echoes_ty1_first_ch(2) - first_echoes_ty1_first_ch(1);
    for k=first_echoes_ty1_first_ch(1):length(shot)
        shot(k)=ceil((k-first_echoes_ty1_first_ch(1) + 1) / k_line_per_shot);
    end
    max(shot)

    %Sort Shot Label in the same way as soring data
    MR_DPsti_shotLabel = MR_DPstiTSE.Copy;
    Data_temp = MR_DPsti_shotLabel.Data;
    [kx_dim, kx_lines] = size(Data_temp);
    shot_label = shot(readdata_ix);
    shot_label_data = repmat(shot_label', [kx_dim, 1]);
    % size(shot_label_data) should have the same size as MR_DPstiTSE.Data
    MR_DPsti_shotLabel.Data = single(shot_label_data);
    MR_DPsti_shotLabel.SortData;
    Ima_Sorted_shot_lable = squeeze(double(MR_DPsti_shotLabel.Data));
    save(current_mat_file,'Ima_Sorted_shot_lable','-append');
    %============Display shot labels=============
    clear shot_label_temp; shot_label_temp(:,:,1,:) = squeeze(Ima_Sorted_shot_lable(32,:,:,2,:));
    figure(605); montage(shot_label_temp,'displayrange',[]); colormap jet; title('shot labels');colorbar
    figure(606); imshow(squeeze(Ima_Sorted_shot_lable(32,:,:,1,1)),[]); colormap jet; colorbar

%-------------------DONE---------------------------

%==============Continue with Image Kspa data extracting
MR_DPstiTSE.SortData;
% >>>>>>>exported kspace data from here
ima_k_spa_data = squeeze(double(MR_DPstiTSE.Data));

MR_DPstiTSE.GridData;
MR_DPstiTSE.RingingFilter;
MR_DPstiTSE.ZeroFill;
MR_DPstiTSE.K2IM;
MR_DPstiTSE.EPIPhaseCorrection;
MR_DPstiTSE.K2IP;
MR_DPstiTSE.GridderNormalization;
MR_DPstiTSE.SENSEUnfold;  %no sense here
MR_DPstiTSE.ConcomitantFieldCorrection;
MR_DPstiTSE.DivideFlowSegments;

MR_DPstiTSE.ZeroFill;
MR_DPstiTSE.FlowPhaseCorrection;
MR_DPstiTSE.RotateImage;
MR_DPstiTSE.ShowData;
ima_data = squeeze(double(MR_DPstiTSE.Data));



Ima_Kx_range(1,:) = MR_DPstiTSE.Parameter.Encoding.KxRange(1,:); %1row for mix = 0, i.e. image data
Ima_Ky_range(1,:) = MR_DPstiTSE.Parameter.Encoding.KyRange(1,:);
Ima_Kz_range(1,:) = MR_DPstiTSE.Parameter.Encoding.KzRange(1,:);


save(current_mat_file,'ima_k_spa_data','Ima_Sorted_shot_lable','ima_data','Ima_Kx_range','Ima_Ky_range','Ima_Kz_range','-append');

%% --------------Navigator processing--------------------
clc; close all;
load(current_mat_file);

nav_kspa_recon_nufft = bart('fft 3',nav_im_recon_nufft); 
[kx, ky, ch, shot] = size(nav_kspa_recon_nufft);
save(current_mat_file,'nav_kspa_recon_nufft','-append');

easy_rigidMotion_parameter_calculation(current_mat_file);   

%--------------END--------------------------------------

%% ------Data Correction Initialization-----------------------------------
clc; close all
load(current_mat_file); 

if(prod(double((size(Ima_Sorted_shot_lable) == size(ima_k_spa_data)))) ~= 1)
    error('k space data matrix has different size as shot label!')
end
% Ideal trajectory
clear traj_matrix_x traj_matrix_y traj_matrix_z

[x_dim, y_dim, z_dim, ch_dim, nsa_dim, other_dim ] = size(ima_k_spa_data);

%find the location of the kspace signal peak
max_kspa_index = find(abs(ima_k_spa_data) == max(abs(ima_k_spa_data(:))));
[kx_center, ky_center, kz_center, ch_center] = ind2sub(size(ima_k_spa_data),max_kspa_index);

for Kz = 1:z_dim
    for Ky = 1:y_dim
        for Kx = 1:x_dim
            traj_matrix_x (Kx,Ky,Kz) = round(Kx+Ima_Kx_range(1)/2-1);  %Kx-kx_center;
            traj_matrix_y (Kx,Ky,Kz) = Ky+Ima_Ky_range(1)-1;  %Ky-ky_center;
            traj_matrix_z (Kx,Ky,Kz) = Kz+Ima_Kz_range(1)-1;  %Kz-kz_center;
        end
    end
end
% [traj_matrix_x, traj_matrix_y, traj_matrix_z] = meshgrid(Ima_Kx_range(1)/2:floor(Ima_Kx_range(2)/2), Ima_Ky_range(1):Ima_Ky_range(2),Ima_Kz_range(1):Ima_Kz_range(2)); 
figure(12);
subplot(131); imagesc(array2mosaic(traj_matrix_x)); axis off
subplot(132); imagesc(array2mosaic(traj_matrix_y)); axis off
subplot(133); imagesc(array2mosaic(traj_matrix_z));colormap jet; axis off



%% correction
% linear_phase_xy(2,ch,shot) is used for trajctory correction in x and y(linear phase error)
% global_phase(ch,shot) is used for global phase error correction
% [nav_kx_dim, nav_ky_dim, nav_ch_dim, nav_shot] = size(nav_kspa_recon_nufft);

nav_shot = max(Ima_Sorted_shot_lable(:));
traj_matrix_x_match = repmat(traj_matrix_x,[1,1,1,nsa_dim]);
traj_matrix_y_match = repmat(traj_matrix_y,[1,1,1,nsa_dim]);
traj_matrix_z_match = repmat(traj_matrix_z,[1,1,1,nsa_dim]);


ch_idx = 6;
traj_matrix_x_cor =traj_matrix_x_match;
traj_matrix_y_cor =traj_matrix_y_match;
traj_matrix_z_cor =traj_matrix_z_match;
DP_Ks_data_all_cor = squeeze(ima_k_spa_data(:,:,:,ch_idx,:));

Sorted_shot_lable_1dyn = squeeze(Ima_Sorted_shot_lable(:,:,:,ch_idx,:));

j = sqrt(-1);
correction_mask = 0 * Sorted_shot_lable_1dyn;

for nr_shot = 1:nav_shot   %dyn for spiral is the same as shot for DP
            
    %  correction points
    indx = find(Sorted_shot_lable_1dyn == nr_shot);  % size(Sorted_shot_lable_1dyn) = kx ky kz
    
%     if((abs(global_phase(ch_idx,nr_shot))>50)||(abs(linear_phase_xy(1,ch_idx,nr_shot))>2*pi)) %global_phase in rad, linear_phase_xy in rad/navigator_pixel
%         disp('discard');
%         DP_Ks_data_all_cor(indx) = 0; %exlcude this shot
    if(0)
    else
        disp('Corrected');
%-------------------------------------------------------------------
%        Perform correction
%-------------------------------------------------------------------
        % linear_phase_xy in rad/navigator_pixel (measured from navigator)
%       traj_matrix_x_cor(indx) = traj_matrix_x_cor(indx) - linear_phase_xy(1,ch_idx,nr_shot) * nav_kx_dim / 2 / pi;  %a pixel in DP kspace is (1/0.152)/m  
%       traj_matrix_y_cor(indx) = traj_matrix_y_cor(indx) - linear_phase_xy(2,ch_idx,nr_shot) * nav_ky_dim / 2 / pi;
%       traj_matrix_z_cor(indx) = traj_matrix_z_cor(indx) - 0;  %no correction for kz for this expierement
%         global phase error in rad
%         DP_Ks_data_all_cor(indx) = DP_Ks_data_all_cor(indx) .* exp(-j*global_phase(ch_idx,nr_shot) / 2); %the simulated introduced phase is doubled by the navigator

% >>>>>>>>>>>>>> use input phase error to correct <<<<<<<<<<<<<<<<<<<<<
        traj_matrix_x_cor(indx) = traj_matrix_x_cor(indx) - linear_phase_error_input_x_in_TSEksp_pixel(nr_shot ); 
        traj_matrix_y_cor(indx) = traj_matrix_y_cor(indx) - linear_phase_error_input_y_in_TSEksp_pixel(nr_shot ); 
%         traj_matrix_z_cor(indx) = traj_matrix_z_cor(indx) - 0;  %no correction for kz for this expierement
%         % global phase error in rad
%         DP_Ks_data_all_cor(indx) = DP_Ks_data_all_cor(indx) .* exp(-j*global_phase_error_input_match_rad(nr_shot) ); 

        

              
        correction_mask(indx) = nr_shot;
        
%         for ix = 1:length(indx)
%             z = ceil(indx(ix)/(labbel_mx_size(1)*labbel_mx_size(2)));
%             
%             mod_xy = mod(indx(ix),(labbel_mx_size(1)*labbel_mx_size(2)));
%             if(mod_xy == 0)
%                mod_xy = (labbel_mx_size(1)*labbel_mx_size(2));
%             end
%             y = ceil(mod_xy/labbel_mx_size(1));
%             x = 1+mod(mod(indx(ix),(labbel_mx_size(1)*labbel_mx_size(2))),labbel_mx_size(1));
%             
%             traj_matrix_x_cor(indx) = traj_matrix_x_cor(x,y,z) - rot(nr_shot,1);
%             traj_matrix_y_cor(x,y,z)= traj_matrix_y_cor(x,y,z) - rot(nr_shot,2);
%             traj_matrix_z_cor(x,y,z) = traj_matrix_z_cor(x,y,z) - rot(nr_shot,3);
% 
%             DP_Ks_data_all_cor(x,y,z) = DP_Ks_data_all_cor(x,y,z) * exp(-j*tran(nr_shot,:));
%             
%             correction_mask(x,y,z) = nr_shot;
%         end
    end
    figure(5)
    subplot(1,2,1);
%     imshow(squeeze(correction_mask(1,:,:,1)),[0 292]);colorbar;colormap jet
    montage(permute(Sorted_shot_lable_1dyn(32,:,:,:),[ 2 3 1 4]),'Displayrange',[0 102]); colormap jet; colorbar;
    subplot(1,2,2);
%     imshow(squeeze(correction_mask(1,:,:,2)),[0 292]);colorbar;colormap jet
    montage(permute(correction_mask(32,:,:,:),[ 2 3 1 4]),'Displayrange',[0 102]); colormap jet; colorbar;
    title(['Shot: ',num2str(nr_shot),'; Datapoints: ',num2str(length(indx))]);
    drawnow();
    pause(0.1);
    
   
    figure(14);
    subplot(331); imagesc(array2mosaic(traj_matrix_x_cor(:,:,:,1))); axis off
    subplot(332); imagesc(array2mosaic(traj_matrix_y_cor(:,:,:,1))); axis off
    subplot(333); imagesc(array2mosaic(traj_matrix_z_cor(:,:,:,1)));colormap jet; axis off
    
%     subplot(334); imagesc(array2mosaic(traj_matrix_x_cor(:,:,:,2))); axis off
%     subplot(335); imagesc(array2mosaic(traj_matrix_y_cor(:,:,:,2))); axis off
%     subplot(336); imagesc(array2mosaic(traj_matrix_z_cor(:,:,:,2)));colormap jet; axis off
% 
%     subplot(337); imagesc(array2mosaic(traj_matrix_x_cor(:,:,:,3))); axis off
%     subplot(338); imagesc(array2mosaic(traj_matrix_y_cor(:,:,:,3))); axis off
%     subplot(339); imagesc(array2mosaic(traj_matrix_z_cor(:,:,:,3)));colormap jet; axis off
    drawnow();
    pause(0.1);
end  
% % only kx in -20:20 ky in -20:20 for recon
idx_outrange = find(traj_matrix_x_cor<-20 | traj_matrix_x_cor > 20);
DP_Ks_data_all_cor(idx_outrange) = 0;
idx_outrange = find(traj_matrix_y_cor<-20 | traj_matrix_y_cor > 20);
DP_Ks_data_all_cor(idx_outrange) = 0;

figure(101);
display_kz = 3;
temp_x = traj_matrix_x_match(:,:,display_kz,:); temp_y = traj_matrix_y_match(:,:,display_kz,:); temp_k = abs(ima_k_spa_data(:,:,display_kz,6,:));
subplot(121); scatter(temp_x(:), temp_y(:),30, temp_k(:),'filled'); title('uncor mag'); xlabel('kx'); ylabel('ky')
temp_x_cor = traj_matrix_x_cor(:,:,display_kz,:); temp_y_cor = traj_matrix_y_cor(:,:,display_kz,:); temp_k_cor = abs(DP_Ks_data_all_cor(:,:,display_kz,:));
subplot(122); scatter(temp_x_cor(:), temp_y_cor(:),30, temp_k_cor(:),'filled');title('cor mag'); xlabel('kx'); ylabel('ky')

figure(102);
temp_x = traj_matrix_x_match(:,:,display_kz,:); temp_y = traj_matrix_y_match(:,:,display_kz,:); temp_k = angle(ima_k_spa_data(:,:,display_kz,6,:));
subplot(121); scatter(temp_x(:), temp_y(:),30, temp_k(:),'filled');title('uncor phase'); xlabel('kx'); ylabel('ky')
temp_x_cor = traj_matrix_x_cor(:,:,display_kz,:); temp_y_cor = traj_matrix_y_cor(:,:,display_kz,:); temp_k_cor = angle(DP_Ks_data_all_cor(:,:,display_kz,:));
subplot(122); scatter(temp_x_cor(:), temp_y_cor(:),30, temp_k_cor(:),'filled');title('cor phase'); xlabel('kx'); ylabel('ky')

figure(103);
subplot(121);
temp_x = traj_matrix_x_match(:,:,display_kz,:); temp_y = traj_matrix_y_match(:,:,display_kz,:); temp_k = Sorted_shot_lable_1dyn(:,:,display_kz);
scatter(temp_x(:), temp_y(:),30, temp_k(:),'filled'); title('shot locations'); xlabel('kx'); ylabel('ky')
subplot(122);
temp_x = traj_matrix_x_cor(:,:,display_kz,:); temp_y = traj_matrix_y_cor(:,:,display_kz,:); temp_k = Sorted_shot_lable_1dyn(:,:,display_kz);
scatter(temp_x(:), temp_y(:),30, temp_k(:),'filled'); title('shot locations'); xlabel('kx'); ylabel('ky')

%reshape to fit bart nufft
traj_matrix_ideal = cat(5, traj_matrix_x_match, traj_matrix_y_match, traj_matrix_z_match);
traj_matrix_ideal = permute(traj_matrix_ideal, [5 1 2 3 4]);
[dim1, dimx, dimy, dimz, dimnas] = size(traj_matrix_ideal);
traj_matrix_ideal_rs = reshape(traj_matrix_ideal, 3, dimx, dimy*dimz*dimnas);

%---------------------------------------------
traj_matrix_cor = cat(5, traj_matrix_x_cor, traj_matrix_y_cor, traj_matrix_z_cor);
traj_matrix_cor = permute(traj_matrix_cor, [5 1 2 3 4]);
traj_matrix_cor_rs = reshape(traj_matrix_cor, 3, dimx, dimy*dimz*dimnas);

DP_Ks_data_all_cor_rs = reshape(DP_Ks_data_all_cor, dimx, dimy*dimz*dimnas);
DP_Ks_data_all_cor_pm = permute(DP_Ks_data_all_cor_rs, [3 1 2]);

DP_Ks_data_all_uncor = squeeze(ima_k_spa_data(:,:,:,ch_idx,:));
DP_Ks_data_all_uncor_rs = reshape(DP_Ks_data_all_uncor, dimx, dimy*dimz*dimnas);
DP_Ks_data_all_uncor_pm = permute(DP_Ks_data_all_uncor_rs, [3 1 2]);

%--------------------------------------------
traj_def = traj_matrix_cor - traj_matrix_ideal;
ttt = squeeze(traj_def(:,32,:,:,1)); figure; montage(permute(ttt, [2 3 4 1]),'Displayrange',[]); colormap jet; colorbar

% =========================================================
%                      BART Recon
% =========================================================
clear DP_uncorrected DP_uncorrected_fftshifted DP_corrected DP_corrected_fftshifted
% Uncorrected
% selected_range = [1:732];
DP_uncorrected(:,:,:,ch_idx) = bart('nufft -i -l0.01', traj_matrix_ideal_rs(:,:,:), DP_Ks_data_all_uncor_pm(:,:,:));
DP_uncorrected_fftshifted = fftshift(fftshift(DP_uncorrected,2),3);
figure(20); montage(permute(abs(DP_uncorrected_fftshifted(:,:,:,ch_idx)),[1 2 4 3]),'displayrange',[]); colormap gray
% Corrected
traj_matrix_cor_rs_shifted = traj_matrix_cor_rs;
% traj_matrix_cor_rs_shifted(2,:,:) = traj_matrix_cor_rs_shifted(2,:,:)+1;
DP_corrected(:,:,:,ch_idx) = bart('nufft -i -l0.01', traj_matrix_cor_rs_shifted(:,:,:), DP_Ks_data_all_cor_pm(:,:,:));
DP_corrected_fftshifted = fftshift(fftshift(DP_corrected,2),3);
figure(21); montage(permute(abs(DP_corrected_fftshifted(:,:,:,ch_idx)),[1 2 4 3]),'displayrange',[]); colormap gray

% PICS
% sens = ones(70, 60, 12);
% DP_correctedPICS(:,:,:,ch_idx) = ((bart('pics -S -r0.001 -t', traj_matrix_cor_rs, DP_Ks_data_all_cor_pm, sens)));
% DP_correctedPICS_fftshifted = fftshift(fftshift(DP_correctedPICS,2),3);
% figure(211); montage(permute(abs(DP_correctedPICS_fftshifted(:,:,:,ch_idx)),[1 2 4 3]),'displayrange',[]); colormap gray



% =========================================================
%                    NUFFT Recon
% =========================================================

% Uncorrected

% reshape trajectory to m*3; data to m*1
traj_nufft_ideal = reshape(traj_matrix_ideal_rs, 3, size(traj_matrix_ideal_rs,2)*size(traj_matrix_ideal_rs,3));
sig_nufft_uncor = col(DP_Ks_data_all_uncor_pm);

% scale trajectory [-pi pi]
scale_factor_x = pi / max(abs(traj_nufft_ideal(1,:)));
scale_factor_y = pi / max(abs(traj_nufft_ideal(2,:)));
scale_factor_z = pi / max(abs(traj_nufft_ideal(3,:)));

traj_nufft_ideal_scaled(:,1) = traj_nufft_ideal(1,:) * scale_factor_x;
traj_nufft_ideal_scaled(:,2) = traj_nufft_ideal(2,:) * scale_factor_y;
traj_nufft_ideal_scaled(:,3) = traj_nufft_ideal(3,:) * scale_factor_z;


clear im_recon_nufft_uncor im_recon_nufft_uncor_fftshifted
x_dim = 64; y_dim = 60; z_dim = 4;
A=nuFTOperator(traj_nufft_ideal_scaled,[x_dim, y_dim, z_dim],ones(x_dim, y_dim, z_dim),6);
% im_recon_nufft_uncor(:,:,:,ch_idx)=regularizedReconstruction(A,sig_nufft_uncor,@L2Norm,0.5,'maxit',25);
im_recon_nufft_uncor(:,:,:,ch_idx)=A'*sig_nufft_uncor;

% im_recon_nufft_uncor_fftshifted = fftshift(fftshift(im_recon_nufft_uncor,2),3);
figure(22); montage(permute(abs(im_recon_nufft_uncor(:,:,:,ch_idx)),[1 2 4 3]),'displayrange',[])

% Corrected
% reshape trajectory to m*3; data to m*1
traj_nufft_cor = reshape(traj_matrix_cor_rs, 3, size(traj_matrix_cor_rs,2)*size(traj_matrix_cor_rs,3));
sig_nufft_cor = col(DP_Ks_data_all_cor_pm);

% scale trajectory [-pi pi]
nufft_dim = round(max(traj_nufft_cor,[],2) - min(traj_nufft_cor,[],2)+1); 


scale_factor_x = pi / max(abs(traj_nufft_cor(1,:)));
scale_factor_y = pi / max(abs(traj_nufft_cor(2,:)));
scale_factor_z = pi / max(abs(traj_nufft_cor(3,:)));

traj_nufft_cor_scaled(:,1) = traj_nufft_cor(1,:) * scale_factor_x;
traj_nufft_cor_scaled(:,2) = traj_nufft_cor(2,:) * scale_factor_y;
traj_nufft_cor_scaled(:,3) = traj_nufft_cor(3,:) * scale_factor_z;



clear im_recon_nufft_cor im_recon_nufft_cor_fftshifted
A=nuFTOperator(traj_nufft_cor_scaled,[nufft_dim(1), nufft_dim(2), nufft_dim(3)],ones(nufft_dim(1), nufft_dim(2), nufft_dim(3)),2);
im_recon_nufft_cor(:,:,:,ch_idx)=regularizedReconstruction(A,sig_nufft_cor,@L2Norm,0.5,'maxit',25);
% im_recon_nufft_cor(:,:,:,ch_idx)=A'*sig_nufft_cor;

% im_recon_nufft_cor_fftshifted = fftshift(fftshift(im_recon_nufft_cor,2),3);
figure(23); montage(permute(abs(im_recon_nufft_cor(:,:,:,ch_idx)),[1 2 4 3]),'displayrange',[])
%% remove empty k-points
empty_idx = find(abs(sig_nufft_cor) ==0 );
sig_nufft_cor_rm_empty = sig_nufft_cor;
traj_nufft_cor_scaled_rm_empty = traj_nufft_cor_scaled;

sig_nufft_cor_rm_empty(empty_idx) = [];
traj_nufft_cor_scaled_rm_empty(empty_idx,:) = [];
A=nuFTOperator(traj_nufft_cor_scaled_rm_empty,[nufft_dim(1), nufft_dim(2), nufft_dim(3)],ones(nufft_dim(1), nufft_dim(2), nufft_dim(3)),2);
im_recon_nufft_cor(:,:,:,ch_idx)=regularizedReconstruction(A,sig_nufft_cor,@L2Norm,0.5,'maxit',25);
figure(24); montage(permute(abs(im_recon_nufft_cor(:,:,:,ch_idx)),[1 2 4 3]),'displayrange',[])


