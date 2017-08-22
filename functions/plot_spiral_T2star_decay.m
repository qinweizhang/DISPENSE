
%%
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_08_20_SND_external_waveform')
%% normal spiral
mr = MRecon('sn_20082017_1428162_12_2_wip_3d_sprialstack_iscanV4.raw');

mr.Parameter.Parameter2Read.typ = 1;
mr.Parameter.Parameter2Read.mix = 0;
mr.ReadData;
k = mr.Data;
size(k)

ch = length(mr.Parameter.Labels.CoilNrs)
AQ_base = mr.Parameter.GetObject('AQ`base'); inter =  AQ_base.GetValue('interval')
figure(4); hold on; plot([1:length(k)]*inter, abs(k(:,15))./max(abs(k(:,15))),'k--')

%% DPsti spiral
disp('spiral Nav. data loading...')
mat_fn ='SND.mat';
fn = 'sn_20082017_1442519_2_2_wip_sc17_dpsti_sosad_linearV4.raw';
 nav_kspa_data_read(fn, mat_fn);

disp('-finished- ');
load(mat_fn)


mr_2 = MRecon(fn);
AQ_DPNav = mr_2.Parameter.GetObject('AQ`DPnav'); inter_2 =  AQ_DPNav.GetValue('interval')

figure(4); hold on; plot([1:length(nav_k_spa_data)].* inter_2, squeeze(abs(nav_k_spa_data(:,3,1))./max(abs(nav_k_spa_data(:,3,1)))),'LineWidth',0.5)