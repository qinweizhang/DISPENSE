
clear; vars_kerry; set(0,'DefaultAxesFontSize',22);

%% Set file names
clear; close all; clc
[f,path]=uigetfile(fullfile(cd,'*.raw'),'M, P,  S','MultiSelect','on');
fn_M = strcat(path,f{1,1});
fn_P = strcat(path,f{1,2});
fn_S = strcat(path,f{1,3});


% =======================================DPnav Trj Measurement M


MR_TSEDPnav_M = MRecon(fn_M);

MR_DPnavspiralM_recon1 = MR_TSEDPnav_M.Copy;
MR_DPnavspiralM_recon1.Parameter.Parameter2Read.typ = 1;
MR_DPnavspiralM_recon1.Parameter.Parameter2Read.mix = 1;  

MR_DPnavspiralM_recon1.ReadData;
MR_DPnavspiralM_recon1.RandomPhaseCorrection;

k_spa_M_data = double(MR_DPnavspiralM_recon1.Data);
[kx, profiles] = size(k_spa_M_data);

ch_nr = length(MR_DPnavspiralM_recon1.Parameter.Labels.CoilNrs);
n_nsa = max(MR_DPnavspiralM_recon1.Parameter.Labels.Index.aver) + 1; 
n_dyn = max(MR_DPnavspiralM_recon1.Parameter.Labels.Index.dyn) + 1;
shots_per_volumn = profiles / ch_nr / n_nsa / n_dyn;

k_spa_M_data = reshape(k_spa_M_data,kx, ch_nr, n_nsa, shots_per_volumn, n_dyn); 
[kx, n_ch,  n_nsa, shots, diffusion_setting] = size(k_spa_M_data)

k_spa_M_data_rm_phase_offset = k_spa_M_data(kx/2+1:end,:,:,:,:); 
k_spa_M1_data_rm_phase_offset = squeeze(k_spa_M_data_rm_phase_offset(:,:,1,:,:));
k_spa_M2_data_rm_phase_offset = squeeze(k_spa_M_data_rm_phase_offset(:,:,2,:,:));
k_spa_M2_data_rm_phase_offset = k_spa_M2_data_rm_phase_offset.* exp(i*pi);

for diffusion_nr = 1:diffusion_setting
    figure(100+diffusion_nr); 
    plot(squeeze(unwrap(angle(k_spa_M1_data_rm_phase_offset(:, 1,:,diffusion_nr))))); title('1 repeatition')
    hold on
    plot(squeeze(unwrap(angle(k_spa_M2_data_rm_phase_offset(:, 1,:,diffusion_nr))))); title('1 repeatition')
    title(['M phase diffusion nr = ',num2str(diffusion_nr)]);
    hold off
    drawnow();
    pause(1);

    figure(104);
    b0_phase = squeeze(unwrap(angle(k_spa_M2_data_rm_phase_offset(:,2,:,1))));
    ec_phase = squeeze(unwrap(angle(k_spa_M2_data_rm_phase_offset(:,2,:,diffusion_nr))));
    hold on;
    plot(ec_phase - b0_phase);  
    drawnow();
    pause(1);
end

k_spa_M1_data_phase_NSA = mean(k_spa_M1_data_rm_phase_offset,3);
k_spa_M2_data_phase_NSA = mean(k_spa_M2_data_rm_phase_offset,3);


%=======================================DPnav Trj Measurement P


MR_TSEDPnav_P = MRecon(fn_P);

MR_DPnavspiralP_recon1 = MR_TSEDPnav_P.Copy;
MR_DPnavspiralP_recon1.Parameter.Parameter2Read.typ = 1;
MR_DPnavspiralP_recon1.Parameter.Parameter2Read.mix = 1;  

MR_DPnavspiralP_recon1.ReadData;
MR_DPnavspiralP_recon1.RandomPhaseCorrection;

k_spa_P_data = double(MR_DPnavspiralP_recon1.Data);
[kx, profiles] = size(k_spa_P_data);

ch_nr = length(MR_DPnavspiralP_recon1.Parameter.Labels.CoilNrs);
n_nsa = max(MR_DPnavspiralP_recon1.Parameter.Labels.Index.aver) + 1; 
n_dyn = max(MR_DPnavspiralP_recon1.Parameter.Labels.Index.dyn) + 1;
shots_per_volumn = profiles / ch_nr / n_nsa / n_dyn;

k_spa_P_data = reshape(k_spa_P_data,kx, ch_nr, n_nsa, shots_per_volumn, n_dyn); 
[kx, n_ch,  n_nsa, shots, diffusion_setting] = size(k_spa_P_data)

k_spa_P_data_rm_phase_offset = k_spa_P_data(kx/2+1:end,:,:,:,:); 
k_spa_P1_data_rm_phase_offset = squeeze(k_spa_P_data_rm_phase_offset(:,:,1,:,:));
k_spa_P2_data_rm_phase_offset = squeeze(k_spa_P_data_rm_phase_offset(:,:,2,:,:));
k_spa_P2_data_rm_phase_offset = k_spa_P2_data_rm_phase_offset.* exp(i*pi);

for diffusion_nr = 1:diffusion_setting
    figure(200++diffusion_nr); 
    plot(squeeze(unwrap(angle(k_spa_P1_data_rm_phase_offset(:, 1,:,diffusion_nr))))); title('1 repeatition')
    hold on
    plot(squeeze(unwrap(angle(k_spa_P2_data_rm_phase_offset(:, 1,:,diffusion_nr))))); title('1 repeatition')
    title(['P phase diffusion nr = ',num2str(diffusion_nr)]);
    hold off
    drawnow();
    pause(1);

    figure(204);
    b0_phase = squeeze(unwrap(angle(k_spa_P2_data_rm_phase_offset(:,2,:,1))));
    ec_phase = squeeze(unwrap(angle(k_spa_P2_data_rm_phase_offset(:,2,:,diffusion_nr))));
    hold on;
    plot(ec_phase - b0_phase);
    drawnow();
    pause(1);
end

k_spa_P1_data_phase_NSA = mean(k_spa_P1_data_rm_phase_offset,3);
k_spa_P2_data_phase_NSA = mean(k_spa_P2_data_rm_phase_offset,3);


%=======================================DPnav Trj Measurement S


MR_TSEDPnav_S = MRecon(fn_S);

MR_DPnavspiralS_recon1 = MR_TSEDPnav_S.Copy;
MR_DPnavspiralS_recon1.Parameter.Parameter2Read.typ = 1;
MR_DPnavspiralS_recon1.Parameter.Parameter2Read.mix = 1;  

MR_DPnavspiralS_recon1.ReadData;
MR_DPnavspiralS_recon1.RandomPhaseCorrection;

k_spa_S_data = double(MR_DPnavspiralS_recon1.Data);
[kx, profiles] = size(k_spa_S_data);

ch_nr = length(MR_DPnavspiralS_recon1.Parameter.Labels.CoilNrs);
n_nsa = max(MR_DPnavspiralS_recon1.Parameter.Labels.Index.aver) + 1; 
n_dyn = max(MR_DPnavspiralS_recon1.Parameter.Labels.Index.dyn) + 1;
shots_per_volumn = profiles / ch_nr / n_nsa / n_dyn;

k_spa_S_data = reshape(k_spa_S_data,kx, ch_nr, n_nsa, shots_per_volumn, n_dyn); 
[kx, n_ch,  n_nsa, shots, diffusion_setting] = size(k_spa_S_data)

k_spa_S_data_rm_phase_offset = k_spa_S_data(kx/2+1:end,:,:,:,:); 
k_spa_S1_data_rm_phase_offset = squeeze(k_spa_S_data_rm_phase_offset(:,:,1,:,:));
k_spa_S2_data_rm_phase_offset = squeeze(k_spa_S_data_rm_phase_offset(:,:,2,:,:));
k_spa_S2_data_rm_phase_offset = k_spa_S2_data_rm_phase_offset.* exp(i*pi);

for diffusion_nr = 1:diffusion_setting
    figure(300++diffusion_nr); 
    plot(squeeze(unwrap(angle(k_spa_S1_data_rm_phase_offset(:, 1,:,diffusion_nr))))); title('1 repeatition')
    hold on
    plot(2*pi+squeeze(unwrap(angle(k_spa_S2_data_rm_phase_offset(:, 1,:,diffusion_nr))))); title('1 repeatition')
    title(['S phase diffusion nr = ',num2str(diffusion_nr)]);
    hold off
    drawnow();
    pause(1);

    figure(304);
    b0_phase = squeeze(unwrap(angle(k_spa_S2_data_rm_phase_offset(:,2,:,1))));
    ec_phase = squeeze(unwrap(angle(k_spa_S2_data_rm_phase_offset(:,2,:,diffusion_nr))));
    hold on;
    plot(ec_phase - b0_phase);
    drawnow();
    pause(1);
end

    k_spa_S1_data_phase_NSA = mean(k_spa_S1_data_rm_phase_offset,3);
    k_spa_S2_data_phase_NSA = mean(k_spa_S2_data_rm_phase_offset,3);

%% %----------Calculation----------------%
im_M1_ksp = squeeze(k_spa_M1_data_phase_NSA);
im_M2_ksp = squeeze(k_spa_M2_data_phase_NSA);
im_P1_ksp = squeeze(k_spa_P1_data_phase_NSA);
im_P2_ksp = squeeze(k_spa_P2_data_phase_NSA);
im_S1_ksp = squeeze(k_spa_S1_data_phase_NSA);
im_S2_ksp = squeeze(k_spa_S2_data_phase_NSA);

 clear im_M1_ksp_phase_unwrap im_M2_ksp_phase_unwrap im_P1_ksp_phase_unwrap im_P2_ksp_phase_unwrap im_S1_ksp_phase_unwrap im_S2_ksp_phase_unwrap
 
sel_ch = [1:ch_nr];
for ch=1:length(sel_ch)
    im_M1_ksp_phase_unwrap(:,ch,:) = unwrap(squeeze(angle(im_M1_ksp(:,sel_ch(ch),:))));
    im_M2_ksp_phase_unwrap(:,ch,:) = unwrap(squeeze(angle(im_M2_ksp(:,sel_ch(ch),:))));
end

sel_ch = [1:ch_nr];
for ch=1:length(sel_ch)
     im_P1_ksp_phase_unwrap(:,ch,:) = unwrap(squeeze(angle(im_P1_ksp(:,sel_ch(ch),:))));
     im_P2_ksp_phase_unwrap(:,ch,:) = unwrap(squeeze(angle(im_P2_ksp(:,sel_ch(ch),:))));
end

sel_ch = [1:ch_nr];
for ch=1:length(sel_ch)
     im_S1_ksp_phase_unwrap(:,ch,:) = unwrap(squeeze(angle(im_S1_ksp(:,sel_ch(ch),:))));
     im_S2_ksp_phase_unwrap(:,ch,:) = unwrap(squeeze(angle(im_S2_ksp(:,sel_ch(ch),:))));
end

diffusion_nr = 2;
b0_diffusion_idx = 1;
mm = im_M1_ksp_phase_unwrap(:,:,diffusion_nr)-im_M2_ksp_phase_unwrap(:,:,b0_diffusion_idx); 
figure(1001); plot(mm(1,:))
pp =im_P1_ksp_phase_unwrap(:,:,diffusion_nr)-im_P2_ksp_phase_unwrap(:,:,b0_diffusion_idx);
figure(1002); plot(pp(1,:))
ss =im_S1_ksp_phase_unwrap(:,:,diffusion_nr)-im_S2_ksp_phase_unwrap(:,:,b0_diffusion_idx);
figure(1003); plot(ss(1,:))

wrap_idx = [1 2 3 4 5 7 8 9 11 ];
im_M1_ksp_phase_unwrap(:,wrap_idx,diffusion_nr) = im_M1_ksp_phase_unwrap(:,wrap_idx,diffusion_nr)-2* pi;
wrap_idx_p = [2 3 4 5  7 8 9 10 11];
im_P1_ksp_phase_unwrap(:,wrap_idx_p,diffusion_nr) = im_P1_ksp_phase_unwrap(:,wrap_idx_p,diffusion_nr) - 2* pi;
wrap_idx_z = [1 2 3 4 5  7 8 9   13];
im_S1_ksp_phase_unwrap(:,wrap_idx_z,diffusion_nr) = im_S1_ksp_phase_unwrap(:,wrap_idx_z,diffusion_nr) - 2* pi;


diffusion_nr = 2;
figure(5);
subplot(311); plot(im_M1_ksp_phase_unwrap(:,:,diffusion_nr)); title('phase M step1: spiral on')
subplot(312); plot(im_M2_ksp_phase_unwrap(:,:,diffusion_nr)); title('phase M step2: spiral off')
subplot(313); plot(im_M1_ksp_phase_unwrap(:,:,diffusion_nr)-im_M2_ksp_phase_unwrap(:,:,b0_diffusion_idx)); title('phase M step1-step2')
figure(6);
subplot(311); plot(im_P1_ksp_phase_unwrap(:,:,diffusion_nr)); title('phase P step1: spiral on')
subplot(312); plot(im_P2_ksp_phase_unwrap(:,:,diffusion_nr)); title('phase P step2: spiral off')
subplot(313); plot(im_P1_ksp_phase_unwrap(:,:,diffusion_nr)-im_P2_ksp_phase_unwrap(:,:,b0_diffusion_idx)); title('phase P step1-step2')
figure(7);
subplot(311); plot(im_S1_ksp_phase_unwrap(:,:,diffusion_nr)); title('phase S step1: spiral on')
subplot(312); plot(im_S2_ksp_phase_unwrap(:,:,diffusion_nr)); title('phase S step2: spiral off')
subplot(313); plot(im_S1_ksp_phase_unwrap(:,:,diffusion_nr)-im_S2_ksp_phase_unwrap(:,:,b0_diffusion_idx)); title('phase S step1-step2')

for diffusion_nr = 1:diffusion_setting
    phi_M(:,diffusion_nr) = squeeze(mean(im_M1_ksp_phase_unwrap(:,:,diffusion_nr)-im_M2_ksp_phase_unwrap(:,:,b0_diffusion_idx),2));
    phi_P(:,diffusion_nr) = squeeze(mean(im_P1_ksp_phase_unwrap(:,:,diffusion_nr)-im_P2_ksp_phase_unwrap(:,:,b0_diffusion_idx),2));
    phi_S(:,diffusion_nr) = squeeze(mean(im_S1_ksp_phase_unwrap(:,:,diffusion_nr)-im_S2_ksp_phase_unwrap(:,:,b0_diffusion_idx),2));
end

d = 15; % in mm
trj_meas_kx = phi_M / d / (2 * pi);
trj_meas_ky = phi_P / d / (2 * pi);
trj_meas_kz = phi_S / d / (2 * pi);


figure(8); plot(trj_meas_kz);

diffusion_nr =1;
figure(7);
tt = length(trj_meas_kx)
plot_range = [1:tt];
subplot(121); plot(trj_meas_kx(plot_range,diffusion_nr)); hold on; plot(trj_meas_ky(plot_range,diffusion_nr),'r');  plot(trj_meas_kz(plot_range,diffusion_nr),'k'); legend('measured kx (mm^-^1)','measured ky (mm^-^1)','measured kz (mm^-^1)');
subplot(122); plot3(trj_meas_kx(plot_range,diffusion_nr), trj_meas_ky(plot_range,diffusion_nr), trj_meas_kz(plot_range,diffusion_nr)); xlabel('measured kx (mm^-^1)'); ylabel('measured ky (mm^-^1)'); zlabel('measured kz (mm^-^1)')

weight = sqrt(diff(squeeze(trj_meas_kx(:,diffusion_nr))).^2 + diff(squeeze(trj_meas_ky(:,diffusion_nr))).^2 + diff(squeeze(trj_meas_kz(:,diffusion_nr))).^2);


save traj_Sc4_5_8.mat trj_meas_kx trj_meas_ky trj_meas_kz
%% Recon with measured Trajecotry 
clear
fn = 'dp_16072017_1706523_13_2_wip_dpsti_sad_2dV4.raw';
MR_TSEDPnav_data = MRecon(fn);

MR_TSEDPnav_data_recon1 = MR_TSEDPnav_data.Copy;

MR_TSEDPnav_data_recon1.Parameter.Parameter2Read.typ = 1;
MR_TSEDPnav_data_recon1.Parameter.Parameter2Read.mix = 1;  %for DPnav Spirals

MR_TSEDPnav_data_recon1.Parameter.Recon.ImmediateAveraging = 'No';
MR_TSEDPnav_data_recon1.ReadData;
MR_TSEDPnav_data_recon1.RandomPhaseCorrection;
MR_TSEDPnav_data_recon1.DcOffsetCorrection;
MR_TSEDPnav_data_recon1.MeasPhaseCorrection;

nav_k_spa_data = double(MR_TSEDPnav_data_recon1.Data);
[kx, profiles] = size(nav_k_spa_data);

ch_nr = length(MR_TSEDPnav_data_recon1.Parameter.Labels.CoilNrs);
n_nsa = max(MR_TSEDPnav_data_recon1.Parameter.Labels.Index.aver) + 1; 
n_dyn = max(MR_TSEDPnav_data_recon1.Parameter.Labels.Index.dyn) + 1;
shots_per_volumn = profiles / ch_nr  / n_dyn;

nav_k_spa_data = reshape(nav_k_spa_data,kx, ch_nr, shots_per_volumn, n_dyn); %for normal Spirals

[kx, n_ch, shots, diffusion_setting] = size(nav_k_spa_data)

nav_k_spa_data = nav_k_spa_data(kx/2+1:end,:,:,:,:); %for DPnav when no pi jump happens
%% DEFAULT RECON 
MR_TSEDPnav_data_recon1.SortData;
MR_TSEDPnav_data_recon1.GridData;
MR_TSEDPnav_data_recon1.RingingFilter;
MR_TSEDPnav_data_recon1.ZeroFill;
MR_TSEDPnav_data_recon1.K2I;

MR_TSEDPnav_data_recon1.GridderNormalization;
MR_TSEDPnav_data_recon1.SENSEUnfold;
% MR_TSEDPnav_data_recon1.PartialFourier; %Cause OUT OF MEMORY
MR_TSEDPnav_data_recon1.ConcomitantFieldCorrection;
MR_TSEDPnav_data_recon1.DivideFlowSegments;
% MR_TSEDPnav_data_recon1.Parameter.Recon.CoilCombination = 'pc'; MR_TSEDPnav_data_recon1.CombineCoils;
% MR_TSE_recon1.GeometryCorrection; %don't do this here
% MR_TSEDPnav_data_recon1.RemoveOversampling; %weird SENSE factor....
MR_TSEDPnav_data_recon1.ZeroFill;
MR_TSEDPnav_data_recon1.FlowPhaseCorrection;
MR_TSEDPnav_data_recon1.RotateImage;

im_data = squeeze(MR_TSEDPnav_data_recon1.Data);
figure(41); montage(abs((im_data(32:96,:,3,:))),'displayrange',[],'size',[4 4]); title('Default recon mag'); 
figure(42); montage(angle((im_data(32:96,:,3,:))),'displayrange',[-pi pi],'size',[4 4]); title('Default recon mag'); colormap jet

MR_TSEDPnav_data_recon1.Parameter.Recon.CoilCombination = 'pc'; MR_TSEDPnav_data_recon1.CombineCoils;
im_data_cc = double(MR_TSEDPnav_data_recon1.Data);
figure(43); imshow(abs(((im_data_cc(32:96,:,3)))),[])
figure(44); imshow(angle(((im_data_cc(32:96,:,3)))),[-pi pi]); colormap jet; colorbar

%% SAVE DATA 
save_fn = 'data_Sc09_3D.mat';
if(exist(save_fn)>0)
    save(save_fn, 'nav_k_spa_data','-append');
else 
    save(save_fn, 'nav_k_spa_data');
end


%% -----BART recon -------%
% =======  LOAD DATA
clear; close all;  load('traj_Sc4_5_6.mat'); load('data_Sc03_3D.mat')


%scle trajectory to match kx_range, ky_range (-8 to 8)
close all;
ignor_kz = 1;  % 1for yes. 0 for no

% figure(801); plot(trj_meas_kz,'o');
% remve_point = [1553:1656 3141:3185 4722:4794 6302:6351 7882:7980 9462:9534 11040:11120 12620:12680];
% trj_meas_kx(remve_point,:) = [];
% trj_meas_ky(remve_point,:) = [];
% trj_meas_kz(remve_point,:) = [];
% k_spa_data(remve_point,:,:,:,:) = [];

diffusion_nr = 1;


skip_point =0;
end_point = 1000;% size(trj_meas_kx,1);
selected_point = [skip_point+1:end_point];
trj_meas_kx_t = trj_meas_kx(selected_point,diffusion_nr);
trj_meas_ky_t = trj_meas_ky(selected_point,diffusion_nr);
trj_meas_kz_t = trj_meas_kz(selected_point,diffusion_nr);

clear trj_meas_kx trj_meas_ky trj_meas_kz
trj_meas_kx = trj_meas_kx_t;
trj_meas_ky = trj_meas_ky_t;
trj_meas_kz = trj_meas_kz_t;

%translate trajectory to isocenter (always starts from (0,0))
% trj_meas_kx = trj_meas_kx-trj_meas_kx(1);
% trj_meas_ky = trj_meas_ky-trj_meas_ky(1);

scale_foctor_kx_ky = max([12/max(abs(trj_meas_kx)),12/max(abs(trj_meas_ky))]);
scale_foctor_kz = 7.5/max(abs(trj_meas_kz));
% scale_foctor_kx_ky = 18; 
trj_meas_kx_scaled = trj_meas_kx * scale_foctor_kx_ky;
trj_meas_ky_scaled = trj_meas_ky * scale_foctor_kx_ky;
trj_meas_kz_scaled = trj_meas_kz * scale_foctor_kz;

% trj_meas_ky_scaled = trj_meas_ky_scaled + 0.61;
% trj_meas_kx_scaled = trj_meas_kx_scaled + 0.255;
figure(8);
subplot(121); plot(trj_meas_kx_scaled); hold on; plot(trj_meas_ky_scaled,'r');plot(trj_meas_kz_scaled,'k'); legend('measured kx (a.u.)','measured ky (a.u.)','measured kz (a.u.)');
subplot(122); plot3(trj_meas_kx_scaled, trj_meas_ky_scaled, trj_meas_kz_scaled); xlabel('measured kx (a.u.)'); ylabel('measured ky (a.u.)'); zlabel('measured kz (a.u.)')

 clear igrid sig_bart
 shot_nr = 1;
 nsa_nr = 1;
 
 if(ignor_kz)
    trj_bart = double(cat(2, trj_meas_kx_scaled, trj_meas_ky_scaled, zeros(size(trj_meas_kx_scaled))));
 else
    trj_bart = double(cat(2, trj_meas_kx_scaled, trj_meas_ky_scaled, trj_meas_kz_scaled));
 end

sig_bart(1,:,1,:) = nav_k_spa_data(selected_point,:,nsa_nr,shot_nr,diffusion_nr);


igrid=((bart('nufft -i -l0.01',trj_bart',sig_bart))); % nufft take trj  in [3,kx_range,[ky_range],[kz_range]]; kspa_data in [1,kx_range,[ky_range],[kz_range],ch]

igrid_rss = bart('rss 8',igrid);
slice_id = 1;
figure(21); montage(abs(igrid(:,:,slice_id,:)),'Displayrange',[],'size',[4 4])
figure(22); montage(angle(igrid(:,:,slice_id,:)),'Displayrange',[-pi pi],'size',[4 4]); colormap 'jet'
figure(23); montage(abs(permute(igrid_rss,[1 2 4 3])),'Displayrange',[])

%% %-----BART PICS-------%

% reconstruct low-resolution image and transform back to k-space

lowres_img = bart('nufft -i -d24:24:4 -l0.01', trj_bart',sig_bart );
lowres_ksp = bart('fft -u 7', lowres_img);

figure(202); montage(abs(lowres_img(:,:,2,:)),'displayrange',[]);

% zeropad to full size
clear fullres_ksp
ksp_zerop = bart('resize -c 0 34 1 32 2 4', lowres_ksp);
fullres_ksp(:,:,:,:) = bart('fft -u 7', igrid);

% ESPIRiT calibration
clear sens
% sens = bart('ecalib -m1', ksp_zerop); %low res sense
sens = bart('ecalib -m1', fullres_ksp); %high res sense

figure(203); montage(abs(sens),'displayrange',[]);

% non-Cartesian parallel imging
reco2 = ((bart('pics -S -r0.001 -t', trj_bart', sig_bart, sens)));

figure(204); imshow(abs(((reco2))),[])
figure(205); imshow(angle(((reco2))),[-pi pi]); colormap jet

save data_DWI_fruit.mat nav_k_spa_data reco2 igrid igrid_rss
%% NUFFT recon. 2D
clear; close all; clc;
fn = 'data_Sc8_2D.mat';
 load('traj_Sc3_4_7.mat'); load(fn);

%scle trajectory to match kx_range, ky_range (-pi to pi)

shot_nr = 1; 
diffusion_nr = 1;
nav_k_spa_data = squeeze(nav_k_spa_data);
skip_point = 0;
end_point = length(trj_meas_kx);
end_point = 2000;
selected_point = [skip_point+1:end_point];

trj_meas_kx_t = squeeze(trj_meas_kx(selected_point,1));
trj_meas_ky_t = squeeze(trj_meas_ky(selected_point,1));
trj_meas_kz_t = squeeze(trj_meas_kz(selected_point,1));

clear trj_meas_kx trj_meas_ky trj_meas_kz
trj_meas_kx = trj_meas_kx_t;
trj_meas_ky = trj_meas_ky_t;
trj_meas_kz = trj_meas_kz_t;

scale_foctor = max(pi/max(abs(trj_meas_kx(:,diffusion_nr))),pi/max(abs(trj_meas_ky(:,diffusion_nr))));
trj_meas_kx_scaled = trj_meas_kx(:,diffusion_nr) * scale_foctor;
trj_meas_ky_scaled = trj_meas_ky(:,diffusion_nr) * scale_foctor;

figure(8);
subplot(121); plot(trj_meas_kx_scaled); hold on; plot(trj_meas_ky_scaled,'r'); legend('measured kx (a.u.)','measured ky (a.u.)');
subplot(122); plot(trj_meas_kx_scaled, trj_meas_ky_scaled); xlabel('measured kx (a.u.)'); ylabel('measured ky (a.u.)')



trj_nufft = double(cat(2, trj_meas_kx_scaled, trj_meas_ky_scaled));

x_recon_dim = 24;
y_recon_dim = 24;
clear im_recon_nufft
diffusion_nr = 1; % recon this diffusion volumn
for shot_nr = 1: size(nav_k_spa_data,3)
% for shot_nr = 1: 3
    sig_kspa = nav_k_spa_data(selected_point,:,shot_nr,diffusion_nr);
    disp(['shot nr: ',num2str(shot_nr)]);
    for ch =1:size(nav_k_spa_data,2)
        sig_nufft = double(sig_kspa(:, ch));
        sig_nufft = sig_nufft';
        A=nuFTOperator(trj_nufft,[x_recon_dim, y_recon_dim],ones(x_recon_dim, y_recon_dim),6);

        % simple inverse
    %     im_recon_nufft(:,:,ch) = A'*sig_nufft';



         %call CG-SENSE with L2-norm regularization
         nav_im_recon_nufft(:,:,ch,shot_nr)=regularizedReconstruction(A,sig_nufft',@L2Norm,0.5,'maxit',25);
    %     im_recon_nufft(:,:,ch)=regularizedReconstruction(A,sig_nufft','maxit',25);
    end
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
%% NUFFT recon.  3D
clear; close all; clc;
load('traj_Sc4_5_6.mat'); load('data_Sc03_3D.mat')
%scle trajectory to match kx_range, ky_range (-pi to pi)
ignore_kz = 0; 
recon_dim  = [24 24 15];

skip_point = 0;
end_point = length(trj_meas_kx);
selected_point = [skip_point+1:end_point];

trj_meas_kx_t = squeeze(trj_meas_kx(selected_point,1));
trj_meas_ky_t = squeeze(trj_meas_ky(selected_point,1));
trj_meas_kz_t = squeeze(trj_meas_kz(selected_point,1));

clear trj_meas_kx trj_meas_ky trj_meas_kz k_spa_data
trj_meas_kx = trj_meas_kx_t';
trj_meas_ky = trj_meas_ky_t';
trj_meas_kz = trj_meas_kz_t';

clear trj_meas_kx_t trj_meas_ky_t trj_meas_kz_t sig_bart_t

scale_foctor_xy = max(pi/max(abs(trj_meas_kx)),pi/max(abs(trj_meas_ky)));
trj_meas_kx_scaled = trj_meas_kx * scale_foctor_xy;
trj_meas_ky_scaled = trj_meas_ky * scale_foctor_xy;

scale_foctor_z = pi/max(abs(trj_meas_kz));

if(ignore_kz == 1)
    scale_foctor_z = 0;
end
trj_meas_kz_scaled = trj_meas_kz * scale_foctor_z;


figure(8);
plot(trj_meas_kx_scaled); hold on; plot(trj_meas_ky_scaled,'r');  plot(trj_meas_kz_scaled,'k');legend('measured kx (a.u.)','measured ky (a.u.)','measured kz (a.u.)');

figure(801);
plot3(trj_meas_kx_scaled, trj_meas_ky_scaled, trj_meas_kz_scaled); title('-pi to pi')


trj_nufft = double(cat(1, trj_meas_kx_scaled, trj_meas_ky_scaled, trj_meas_kz_scaled));
trj_nufft = trj_nufft';
clear im_recon_nufft nav_im_recon_nufft

diffusion_nr =2;

% for shot_nr = 1: size(nav_k_spa_data,3)
for shot_nr = 1:1
    
    shot_nr
    
    sig_kspa = nav_k_spa_data(selected_point,:,shot_nr,diffusion_nr);
    
    for ch =1:size(nav_k_spa_data,2)
        sig_nufft = double(sig_kspa(:, ch));
        sig_nufft = sig_nufft';
        A=nuFTOperator(trj_nufft,recon_dim,ones(recon_dim),6);
        
        % simple inverse
        %     nav_im_recon_nufft(:,:,:,ch) = A'*sig_nufft';
        
        
        
        %call CG-SENSE with L2-norm regularization
        %      nav_im_recon_nufft(:,:,:,ch)=regularizedReconstruction(A,sig_nufft',@L2Norm,0.5,'maxit',10);
        nav_im_recon_nufft(:,:,:,ch,shot_nr)=regularizedReconstruction(A,sig_nufft','maxit',2);
    end
    % im_recon_nufft = flipdim(flipdim(im_recon_nufft,1),2);
end

display_shot_nr = 1;
figure(35); montage(permute(abs(squeeze(nav_im_recon_nufft(:,:,:,6,display_shot_nr))), [1 2 4 3]),'displayrange',[]); title('all slices ch6')
figure(36); montage(-1* permute(angle(squeeze(nav_im_recon_nufft(:,:,:,6,display_shot_nr))), [1 2 4 3]),'displayrange',[-pi pi]); colormap jet; title('all slices ch6')

nav_im_recon_nufft_rss = squeeze(bart('rss 8', nav_im_recon_nufft));
figure(37); montage(permute(abs(nav_im_recon_nufft_rss), [1 2 4 3]),'displayrange',[]); title('all slices rss')
figure(38); montage(permute(abs(nav_im_recon_nufft_rss), [2 3 4 1]),'displayrange',[]); title('all slices rss')
