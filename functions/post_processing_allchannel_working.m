        %% set Matlab path and TOOLBOX_PATH environment variable
clear; close all;
clc
current_mat_file = 'data_Sc6_2D.mat';
raw_data_fn = 'dp_10052017_1633195_6_2_wip3ddpnavlinearexperiment1senseV4.raw';
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


MR_DPstiTSE_raw = MRecon(raw_data_fn);

MR_DPstiTSE = MR_DPstiTSE_raw.Copy;

MR_DPstiTSE.Parameter.Parameter2Read.typ = 1;
MR_DPstiTSE.Parameter.Parameter2Read.mix = 0;  
MR_DPstiTSE.Parameter.Recon.ImmediateAveraging = 'No';
MR_DPstiTSE.ReadData;
MR_DPstiTSE.RandomPhaseCorrection;
% MR_DPstiTSE.RemoveOversampling;
MR_DPstiTSE.PDACorrection;
% MR_DPstiTSE.DcOffsetCorrection;
MR_DPstiTSE.MeasPhaseCorrection;

%--------------creat a label for shot number--------------&
    dyn=MR_DPstiTSE.Parameter.Labels.Index.dyn;
    typ=MR_DPstiTSE.Parameter.Labels.Index.typ;
    mix=MR_DPstiTSE.Parameter.Labels.Index.mix; 
    nsa=MR_DPstiTSE.Parameter.Labels.Index.aver;
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

    %-------------match labels ky kz shot to data
    ky_all = MR_DPstiTSE.Parameter.Labels.Index.ky;
    kz_all = MR_DPstiTSE.Parameter.Labels.Index.kz;
    
    ky_matched = ky_all(readdata_ix);
    kz_matched = kz_all(readdata_ix);
    shot_matched = shot(readdata_ix);
    nsa_matched = nsa(readdata_ix);
   

    %============Display shot labels=============
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
    clear shot_label_temp; shot_label_temp(:,:,1,:) = squeeze(Ima_Sorted_shot_lable(32,:,:,2,:));
    figure(605); montage(shot_label_temp,'displayrange',[]); colormap jet; title('shot labels');colorbar
    figure(606); imshow(squeeze(Ima_Sorted_shot_lable(32,:,:,1,1)),[]); colormap jet; colorbar
   
%-------------------DONE---------------------------

% >>>>>>>exported kspace data from here
% ima_k_spa_data_temp = squeeze(double(MR_DPstiTSE.Data));
ima_k_spa_data = squeeze(double(MR_DPstiTSE.Data));
% for prof_nr = 1:length(ky_all)
%     ima_k_spa_data = ima_k_spa_data_temp(:,prof_nr).* (-1)^(ky_all(prof_nr) + kz_all(prof_nr));
% end
 


%==============Continue with Image Kspa data extracting
MR_DPstiTSE.SortData;
kspa_sorted = squeeze(double(MR_DPstiTSE.Data));
MR_DPstiTSE.GridData;
MR_DPstiTSE.RingingFilter;
MR_DPstiTSE.ZeroFill;
MR_DPstiTSE.K2IM;
MR_DPstiTSE.EPIPhaseCorrection;
MR_DPstiTSE.K2IP;
MR_DPstiTSE.GridderNormalization;

% % % --------------Calculate SENSE object-------------
% sense_ref = 'dp_25072017_1728337_1000_7_wipsenserefscanV4.raw';
% MR2 = raw_data_fn;
% coil_survey = 'dp_25072017_1725335_1000_2_wipcoilsurveyscanV4.raw';
% MR_sense = MRsense(sense_ref, MR2, coil_survey);
% MR_sense.Mask = 1;
% MR_sense.MatchTargetSize = 1;
% MR_sense.Perform;
% % % ----------------------end-----------------------
% MR_DPstiTSE.Parameter.Recon.Sensitivities = MR_sense;

MR_DPstiTSE.SENSEUnfold;  %no sense here
MR_DPstiTSE.ConcomitantFieldCorrection;
MR_DPstiTSE.DivideFlowSegments;

MR_DPstiTSE.ZeroFill;
MR_DPstiTSE.FlowPhaseCorrection;
MR_DPstiTSE.RotateImage;
MR_DPstiTSE.ShowData;
ima_data = squeeze(double(MR_DPstiTSE.Data));

ch_dim = length(MR_DPstiTSE.Parameter.Labels.CoilNrs);
save(current_mat_file,'ima_k_spa_data','ky_matched','kz_matched','shot_matched', 'ch_dim','kspa_sorted','-append');
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

if(length(ima_k_spa_data)~=length(shot_matched))|(length(ima_k_spa_data)~=length(ky_matched))
    error('k space data matrix has different size as shot label!')
end
% Sort channel
kx_dim = size(ima_k_spa_data, 1);
profiles = size(ima_k_spa_data, 2)/ch_dim;
ima_k_spa_data_sort_ch = reshape(ima_k_spa_data, kx_dim, ch_dim,  profiles);
ky_matched_profiles = ky_matched(1:ch_dim:end);
kz_matched_profiles = kz_matched(1:ch_dim:end);
kx_matched_profiles = repmat([-1* kx_dim/2 :  kx_dim/2-1],length(ky_matched_profiles),1);
shot_matched_profiles = shot_matched(1:ch_dim:end);

%% correction
% linear_phase_xy(2,ch,shot) is used for trajctory correction in x and y(linear phase error)
% global_phase(ch,shot) is used for global phase error correction

nav_shot = max(shot_matched_profiles(:));

traj_1D_x_cor =kx_matched_profiles; % [klines, kx_points]
traj_1D_y_cor =ky_matched_profiles; % [klines, 1]
traj_1D_z_cor =kz_matched_profiles; % [klines, 1]

 %pi shift of raw data? YES! PLEASE
for prof = 1:size(ima_k_spa_data_sort_ch, 3)
    ima_k_spa_data_sort_ch_pi_shifted(:,:,prof) = ima_k_spa_data_sort_ch(:,:,prof) .* ...
        (-1)^(double(ky_matched_profiles(prof))+double(kz_matched_profiles(prof)) );

end

DP_Ks_data_all_cor = ima_k_spa_data_sort_ch_pi_shifted; % [kx_points, ch, klines]


j = sqrt(-1);
dummyshot = 0; 
for nr_shot = 1:nav_shot   %dyn for spiral is the same as shot for DP
            
    nr_shot
    %  correction points
    indx = find(shot_matched_profiles == nr_shot);  % size(Sorted_shot_lable_1dyn_all_channel) = kx ky kz
    
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
%       traj_1D_x_cor(indx,:) = traj_1D_x_cor(indx,:) - linear_phase_xy(1,1,nr_shot) * nav_kx_dim / 2 / pi;  %a pixel in DP kspace is (1/0.152)/m  
%       traj_1D_y_cor(indx) = traj_1D_y_cor(indx) - linear_phase_xy(2,1,nr_shot) * nav_ky_dim / 2 / pi;
%       traj_1D_z_cor(indx) = traj_1D_z_cor(indx) - 0;  %no correction for kz for this expierement
%         global phase error in rad
%         DP_Ks_data_all_cor(:,:,indx) = DP_Ks_data_all_cor(:,:,indx) .* exp(-j*global_phase(1,nr_shot) ); %the simulated introduced phase is doubled by the navigator

% >>>>>>>>>>>>>> use input phase error to correct <<<<<<<<<<<<<<<<<<<<<
%         traj_1D_x_cor(indx,:) = traj_1D_x_cor(indx,:) - linear_phase_error_input_x_in_TSEksp_pixel(dummyshot + nr_shot ) ; 
%         traj_1D_y_cor(indx) = traj_1D_y_cor(indx) - linear_phase_error_input_y_in_TSEksp_pixel(dummyshot + nr_shot ); 
%         traj_1D_z_cor(indx) = traj_1D_z_cor(indx) - linear_phase_error_input_z_in_TSEksp_pixel(dummyshot + nr_shot  ); 
%         % global phase error in rad
        DP_Ks_data_all_cor(:,:,indx) = DP_Ks_data_all_cor(:,:,indx) .* exp(-j*global_phase_error_input_match_rad(dummyshot + nr_shot ) ); 
%         if(nr_shot~=3)
%             DP_Ks_data_all_cor(:,:,indx) = 0;
%         end
            
              
       
    end
    
end  


%% RECON
% =========================================================
%                      BART Recon
% =========================================================
clear DP_uncorrected DP_uncorrected_fftshifted DP_corrected DP_corrected_fftshifted
% Uncorrected
traj_matrix_ideal_rs = cat(3, kx_matched_profiles', repmat(ky_matched_profiles', [ kx_dim 1]), repmat(kz_matched_profiles', [ kx_dim 1]));
traj_matrix_ideal_rs = permute(traj_matrix_ideal_rs,[3 1 2]); % [3 kx profiles]
DP_Ks_data_all_uncor_pm = permute(ima_k_spa_data_sort_ch_pi_shifted, [4 1 3 2]); %[1 kx profiles ch]

DP_uncorrected = bart('nufft -i -l0.01 ', traj_matrix_ideal_rs, DP_Ks_data_all_uncor_pm);
% DP_uncorrected_fftshifted = fftshift(fftshift(DP_uncorrected,2),3);

% figure(20); montage(permute(abs(DP_uncorrected),[1 2 4 3]),'displayrange',[]); colormap gray
ch_idx = 6;
figure(20); montage(permute(abs(DP_uncorrected(:,:,:,ch_idx)),[1 2 4 3]),'displayrange',[]); colormap gray
figure(201); montage(abs(DP_uncorrected(:,:,2,:)),'displayrange',[]); colormap gray

% Corrected
traj_matrix_cor_rs = cat(3, traj_1D_x_cor', repmat(traj_1D_y_cor', [ kx_dim 1]), repmat(traj_1D_z_cor', [ kx_dim 1]));
traj_matrix_cor_rs = permute(traj_matrix_cor_rs,[3 1 2]); % [3 kx profiles]
DP_Ks_data_all_cor_pm = permute(DP_Ks_data_all_cor, [4 1 3 2]); %[1 kx profiles ch]

DP_corrected = bart('nufft -i -l0.1 ', traj_matrix_cor_rs, DP_Ks_data_all_cor_pm);
% DP_corrected_fftshifted = fftshift(fftshift(DP_corrected,2),3);
figure(21); montage(permute(abs(DP_corrected(:,:,:,ch_idx)),[1 2 4 3]),'displayrange',[]); colormap gray
figure(211); montage(abs(DP_corrected(:,:,2,:)),'displayrange',[]); colormap gray

% PICS
%%--------SCANER SENSE-----------------------------------------
%   sens = gen_sense_map( raw_data_fn, [ 68 60 2]);
%   sens_shifted = fftshift(fftshift(sens,2),3);
%%--------ESPiRIT SENSE----------------------------------------
%     [x, y, z, c] = size(DP_corrected_fftshifted);
%     DP_corrected_fftshifted_2FOV = bart('resize -c 0 132 1 116 2 3', DP_corrected_fftshifted); %bigger FOV
%     for slice =1:z
%         
%         fullres_ksp = bart('fft -u 7', DP_corrected_fftshifted_2FOV(:,:,slice,:));
% 
%         sens_ecalib(:,:,slice,:) = bart('ecalib -m1', fullres_ksp); %high res sense
%     end

DP_correctedPICS = ((bart('pics -S  -r0.01 -t', traj_matrix_cor_rs, DP_Ks_data_all_cor_pm, fftshift(fftshift(sens_ecalib,2),3))));
DP_correctedPICS_fftshifted = fftshift(fftshift(DP_correctedPICS,2),3);
figure(211); montage(permute(abs(DP_correctedPICS_fftshifted),[1 2 4 3]),'displayrange',[]); colormap gray


%%

% =========================================================
%                    NUFFT Recon
% =========================================================

shot_matched_profiles_rs = repmat(shot_matched_profiles', [kx_dim, 1]);
shot_all_points_1d = col(shot_matched_profiles_rs);

%columnized trj and data for nufft recon 
% reshape trajectory to m*3; data to m*1

    % uncorrected 
    kx_matched_profiles_rs = col(kx_matched_profiles');
    ky_matched_profiles_rs = col(repmat(ky_matched_profiles', [ kx_dim 1]));
    kz_matched_profiles_rs = col(repmat(kz_matched_profiles', [ kx_dim 1]));
    traj_nufft_uncor = double(cat(1, kx_matched_profiles_rs', ky_matched_profiles_rs', kz_matched_profiles_rs'));
    
    sig_nufft_uncor = reshape(permute(ima_k_spa_data_sort_ch_pi_shifted, [1 3 2]), kx_dim*profiles, ch_dim);

    % corrected 
    traj_1D_x_cor_rs = col(traj_1D_x_cor');
    traj_1D_y_cor_rs = col(repmat(traj_1D_y_cor', [ kx_dim 1]));
    traj_1D_z_cor_rs = col(repmat(traj_1D_z_cor', [ kx_dim 1]));
    traj_nufft_cor = double(cat(1, traj_1D_x_cor_rs', traj_1D_y_cor_rs', traj_1D_z_cor_rs'));
    
    
    sig_nufft_cor = reshape(permute(DP_Ks_data_all_cor, [1 3 2]), kx_dim*profiles, ch_dim);
%>>>>>>choose recon range
% recon_range = find((shot_all_points_1d <= 222)&(mod(shot_all_points_1d, 6)<=2) ); %use the first half
recon_range = find((shot_all_points_1d <= 222) ); %use the first half

% recon_range = find(mod(shot_all_points_1d,10)==1 ); %use the 3rd shot only
% recon_range = [1:length(shot_all_points_1d)]; %all shots

% Uncorrected

nufft_dim = round(max(traj_nufft_uncor,[],2) - min(traj_nufft_uncor,[],2)+1); 

% scale trajectory [-pi pi]
scale_factor_x = pi / max(abs(traj_nufft_uncor(1,:)));
scale_factor_y = pi / max(abs(traj_nufft_uncor(2,:)));
scale_factor_z = pi / max(abs(traj_nufft_uncor(3,:)));

traj_nufft_ideal_scaled(:,1) = traj_nufft_uncor(1,:) * scale_factor_x;
traj_nufft_ideal_scaled(:,2) = traj_nufft_uncor(2,:) * scale_factor_y;
traj_nufft_ideal_scaled(:,3) = traj_nufft_uncor(3,:) * scale_factor_z;
traj_nufft_ideal_scaled = double(traj_nufft_ideal_scaled);

clear im_recon_nufft_uncor im_recon_nufft_uncor_fftshifted
x_dim = nufft_dim(1); y_dim = nufft_dim(2); z_dim = nufft_dim(3);

oversample_factor = 3; 
sense_map = 0; %1 available 0: not available
if(sense_map == 1)  %use sense imformation
    sens_nuff_uncor = gen_sense_map( raw_data_fn, [x_dim, y_dim, z_dim]); 
    sens_nuff_uncor_shifted = fftshift(fftshift(sens_nuff_uncor,2),3);
    A=nuFTOperator(traj_nufft_ideal_scaled(recon_range,:),[x_dim, y_dim, z_dim],sens_nuff_uncor_shifted,oversample_factor);
else
    A=nuFTOperator(traj_nufft_ideal_scaled(recon_range,:),[x_dim, y_dim, z_dim],ones(x_dim, y_dim, z_dim),oversample_factor);
end
% im_recon_nufft_uncor=regularizedReconstruction(A,col(sig_nufft_uncor),@L2Norm,0.5,'maxit',25);
im_recon_nufft_uncor=A'*col(sig_nufft_uncor(recon_range,:));

% im_recon_nufft_uncor_fftshifted = fftshift(fftshift(im_recon_nufft_uncor,2),3);
figure(22); montage(permute(abs(im_recon_nufft_uncor),[1 2 4 3]),'displayrange',[]); title('uncor')

% Corrected



% scale trajectory [-pi pi]
nufft_dim = round(max(traj_nufft_cor,[],2) - min(traj_nufft_cor,[],2)+1); 


scale_factor_x = pi / max(abs(traj_nufft_cor(1,:)));
scale_factor_y = pi / max(abs(traj_nufft_cor(2,:)));
scale_factor_z = pi / max(abs(traj_nufft_cor(3,:)));

traj_nufft_cor_scaled(:,1) = traj_nufft_cor(1,:) * scale_factor_x;
traj_nufft_cor_scaled(:,2) = traj_nufft_cor(2,:) * scale_factor_y;
traj_nufft_cor_scaled(:,3) = traj_nufft_cor(3,:) * scale_factor_z;

clear im_recon_nufft_cor im_recon_nufft_cor_fftshifted

if(sense_map == 1)  %use sense imformation

    sens_nuff_cor = gen_sense_map( raw_data_fn, [nufft_dim(1), 2 * nufft_dim(2), nufft_dim(3)]); 
    sens_nuff_cor_shifted = fftshift(fftshift(sens_nuff_cor,2),3);
    A=nuFTOperator(traj_nufft_cor_scaled(recon_range,:),[nufft_dim(1), 2 * nufft_dim(2), nufft_dim(3)],sens_nuff_cor_shifted,oversample_factor);
else
%     A=nuFTOperator(traj_nufft_cor_scaled,[nufft_dim(1), nufft_dim(2), nufft_dim(3)],ones(nufft_dim(1), nufft_dim(2), nufft_dim(3),ch_dim),oversample_factor);
    A=nuFTOperator(traj_nufft_cor_scaled(recon_range,:),[nufft_dim(1), nufft_dim(2), nufft_dim(3)],ones(nufft_dim(1), nufft_dim(2), nufft_dim(3)),oversample_factor);
end



% im_recon_nufft_cor=regularizedReconstruction(A,col(sig_nufft_cor(recon_range,:)),@L2Norm,0.5,'maxit',25);
im_recon_nufft_cor=A'*col(sig_nufft_cor(recon_range,:));

im_recon_nufft_cor_fftshifted = fftshift(fftshift(im_recon_nufft_cor,2),3);
figure(24); montage(permute(abs(im_recon_nufft_cor ),[1 2 4 3]),'displayrange',[]); title('cor')

kspa_recon_nuff_cor = bart('fft -u 7',im_recon_nufft_cor);
%% DISPLAY K SPACE

disp_ch = 6; disp_kz = 0;
figure(101);

disp_range = intersect(find(traj_nufft_cor(3,:) == disp_kz), recon_range);

temp_x = traj_nufft_uncor(1,disp_range); temp_y = traj_nufft_uncor(2,disp_range); temp_k = abs(sig_nufft_uncor(disp_range, disp_ch));
subplot(121); scatter(temp_x(:), temp_y(:),30, temp_k(:),'filled'); title('uncor mag'); xlabel('kx'); ylabel('ky')

disp_range = intersect(find(traj_nufft_cor(3,:) == disp_kz), recon_range);
temp_x_cor = traj_nufft_cor(1,disp_range); temp_y_cor = traj_nufft_cor(2,disp_range); temp_k_cor = abs(sig_nufft_cor(disp_range, disp_ch));
subplot(122); scatter(temp_x_cor(:), temp_y_cor(:),30, temp_k_cor(:),'filled');title('cor mag'); xlabel('kx'); ylabel('ky')

figure(102);
disp_range = intersect(find(traj_nufft_uncor(3,:) == disp_kz), recon_range);

temp_x = traj_nufft_uncor(1,disp_range); temp_y = traj_nufft_uncor(2,disp_range); temp_k = angle(sig_nufft_uncor(disp_range, disp_ch));
subplot(121); scatter(temp_x(:), temp_y(:),30, temp_k(:),'filled'); title('uncor phase'); xlabel('kx'); ylabel('ky')

disp_range = intersect(find(traj_nufft_cor(3,:) == disp_kz), recon_range);
temp_x_cor = traj_nufft_cor(1,disp_range); temp_y_cor = traj_nufft_cor(2,disp_range); temp_k_cor = angle(sig_nufft_cor(disp_range, disp_ch));
subplot(122); scatter(temp_x_cor(:), temp_y_cor(:),30, temp_k_cor(:),'filled');title('cor phase'); xlabel('kx'); ylabel('ky')


figure(103);
disp_range = intersect(find(traj_nufft_uncor(3,:) == disp_kz), recon_range);


temp_x = traj_nufft_uncor(1,disp_range); temp_y = traj_nufft_uncor(2,disp_range); temp_k = shot_all_points_1d(disp_range);
subplot(121); scatter(temp_x(:), temp_y(:),30, temp_k(:),'filled'); title('uncor shot location'); xlabel('kx'); ylabel('ky')

disp_range = intersect(find(traj_nufft_cor(3,:) == disp_kz), recon_range);
temp_x_cor = traj_nufft_cor(1,disp_range); temp_y_cor = traj_nufft_cor(2,disp_range); temp_k_cor =  shot_all_points_1d(disp_range);
subplot(122); scatter(temp_x_cor(:), temp_y_cor(:),30, temp_k_cor(:),'filled');title('cor shot location'); xlabel('kx'); ylabel('ky')


figure(104);
subplot(121); imshow(squeeze(abs(kspa_recon_nuff_cor(:,:,2)))',[]);title('nufft recon mag'); xlabel('kx'); ylabel('ky'); colormap jet

subplot(122); imshow(squeeze(angle(kspa_recon_nuff_cor(:,:,2)))',[]);title('nufft recon mag'); xlabel('kx'); ylabel('ky'); colormap jet
