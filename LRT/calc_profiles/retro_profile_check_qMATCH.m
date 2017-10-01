%% measured profile
% cd('/home/qzhang/lood_storage/divi/Ima/parrec/Jasper/LRT/Low_Rank_2017_07_28');
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_09_14\2017_09_14\qM_97389')
MR=MRecon('qm_14092017_1745453_2_2_wip_qmatch_te70V4.raw')
coil_nr = length(MR.Parameter.Labels.CoilNrs);
read_idx = find(MR.Parameter.Labels.Index.typ == 1);
ky = MR.Parameter.Labels.Index.ky;
kz = MR.Parameter.Labels.Index.kz;
dyn = MR.Parameter.Labels.Index.dyn;
rf = MR.Parameter.Labels.Index.rf;

aver= MR.Parameter.Labels.Index.aver;
aver= repmat([0:199],8); aver=aver(:);
aver=repmat(aver,[36]);aver=aver(:); 

ky_read = ky(read_idx); ky_read = ky_read(1:coil_nr:end);
kz_read = kz(read_idx); kz_read = kz_read(1:coil_nr:end);
dyn_read = dyn(read_idx); dyn_read = dyn_read(1:coil_nr:end);
rf_read = rf(read_idx); rf_read = rf_read(1:coil_nr:end);
aver_read = aver(read_idx); aver_read = aver_read(1:coil_nr:end);


ky_range = (max(ky_read)- min(ky_read) + 1)
kz_range = (max(kz_read)- min(kz_read) + 1)
dyn_range = (max(dyn_read)- min(dyn_read) + 1)
rf_range = (max(rf_read)- min(rf_read) + 1)
aver_range = (max(aver_read)- min(aver_read) + 1)

kspa_frames = zeros(ky_range,kz_range,1,dyn_range,aver_range);
% kspa_frames = zeros(ky_range,kz_range,rf_range,dyn_range);

for ii = 1:length(ky_read)
    ky_cor = ky_read(ii) - min(ky_read) + 1;
    kz_cor = kz_read(ii) - min(kz_read) + 1;
    dyn_cor = dyn_read(ii) - min(dyn_read) + 1;
    rf_cor = rf_read(ii) - min(rf_read) + 1;
    aver_cor = aver_read(ii) - min(aver_read) + 1;

    kspa_frames(ky_cor, kz_cor, 1, dyn_cor,aver_cor) = 1 + kspa_frames(ky_cor, kz_cor, 1, dyn_cor,aver_cor);
%     kspa_frames(ky_cor, kz_cor, rf_cor, dyn_cor) = 1 + kspa_frames(ky_cor, kz_cor, rf_cor, dyn_cor);
end

figure(1); montage(squeeze(kspa_frames),'displayrange',[0 3],'size',[6 9]); title('actual sampling pattern'); colormap jet;

% kkk = reshape(kspa_frames, 127, 127, 1, 300);
% figure(1); montage(kkk,'displayrange',[0 3],'size',[5 60]); title('actual sampling pattern'); colormap jet;
%% input profile
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/CSENSE_profiles');
profile = dlmread('LRT_TSE_T2prep_64_64_9_6_r0_l1_bCtr5_sCtr2_us0.05.dat');

nr_profiles_per_dyn = coil_nr0;
ky_input_range = (max(profile(:,1))- min(profile(:,1)) + 1)
kz_input_range = (max(profile(:,2))- min(profile(:,2)) + 1)
dyn_input_range = ceil(length(profile)/nr_profiles_per_dyn)

profile_frames = zeros(ky_input_range,kz_input_range,1,dyn_input_range);

for ii = 1:length(profile)
    ky_input_cor = profile(ii,1) - min(profile(:,1)) + 1;
    kz_input_cor = profile(ii,2) - min(profile(:,2)) + 1;
    dyn_input_cor = ceil(ii / nr_profiles_per_dyn);
    
    profile_frames(ky_input_cor, kz_input_cor, 1, dyn_input_cor) = 1 + profile_frames(ky_input_cor, kz_input_cor, 1, dyn_input_cor);
end
figure(2); montage(profile_frames,'displayrange',[0 3],'size',[6 9]); title('input sampling pattern'); colormap jet;
    
    