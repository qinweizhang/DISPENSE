function spira_3D_trjectory_calculation_normal_6steps(trj_save_fn, d)

% d: offcenter trajectory measurement distance (in mm)

[f,path]=uigetfile(fullfile(cd,'*.raw'),'M, P,  S','MultiSelect','on');
fn_M1 = strcat(path,f{1,1});
fn_M2 = strcat(path,f{1,2});
fn_P1 = strcat(path,f{1,3});
fn_P2 = strcat(path,f{1,4});
fn_S1 = strcat(path,f{1,5});
fn_S2 = strcat(path,f{1,6});

% =======================================DPnav Trj Measurement M

MR_DPnavspiralM1 = MRecon(fn_M1);
MR_DPnavspiralM1.Parameter.Parameter2Read.typ = 1;
MR_DPnavspiralM1.Parameter.Parameter2Read.mix = 0;
MR_DPnavspiralM1.ReadData;
MR_DPnavspiralM1.RandomPhaseCorrection;
k_spa_M1_data = double(MR_DPnavspiralM1.Data);
[kx, profiles] = size(k_spa_M1_data);
ch_nr_m = length(MR_DPnavspiralM1.Parameter.Labels.CoilNrs);
n_nsa = max(MR_DPnavspiralM1.Parameter.Labels.Index.aver) + 1;
n_dyn = max(MR_DPnavspiralM1.Parameter.Labels.Index.dyn) + 1;
shots_per_volumn = profiles / ch_nr_m / n_nsa / n_dyn;
k_spa_M1_data = reshape(k_spa_M1_data,kx, ch_nr_m, n_nsa, shots_per_volumn, n_dyn);
[kx, n_ch,  n_nsa, shots, dynamics] = size(k_spa_M1_data)
k_spa_M1_data_rm_phase_offset = k_spa_M1_data;


MR_DPnavspiralM2 = MRecon(fn_M2);
MR_DPnavspiralM2.Parameter.Parameter2Read.typ = 1;
MR_DPnavspiralM2.Parameter.Parameter2Read.mix = 0;
MR_DPnavspiralM2.ReadData;
MR_DPnavspiralM2.RandomPhaseCorrection;
k_spa_M2_data = double(MR_DPnavspiralM2.Data);
[kx, profiles] = size(k_spa_M2_data);
ch_nr_m = length(MR_DPnavspiralM2.Parameter.Labels.CoilNrs);
n_nsa = max(MR_DPnavspiralM2.Parameter.Labels.Index.aver) + 1;
n_dyn = max(MR_DPnavspiralM2.Parameter.Labels.Index.dyn) + 1;
shots_per_volumn = profiles / ch_nr_m / n_nsa / n_dyn;
k_spa_M2_data = reshape(k_spa_M2_data,kx, ch_nr_m, n_nsa, shots_per_volumn, n_dyn);
[kx, n_ch,  n_nsa, shots, dynamics] = size(k_spa_M2_data)
k_spa_M2_data_rm_phase_offset = k_spa_M2_data;


for diffusion_nr = 1:dynamics
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


MR_DPnavspiralP1 = MRecon(fn_P1);
MR_DPnavspiralP1.Parameter.Parameter2Read.typ = 1;
MR_DPnavspiralP1.Parameter.Parameter2Read.mix = 0;
MR_DPnavspiralP1.ReadData;
MR_DPnavspiralP1.RandomPhaseCorrection;
k_spa_P1_data = double(MR_DPnavspiralP1.Data);
[kx, profiles] = size(k_spa_P1_data);
ch_nr_p = length(MR_DPnavspiralP1.Parameter.Labels.CoilNrs);
n_nsa = max(MR_DPnavspiralP1.Parameter.Labels.Index.aver) + 1;
n_dyn = max(MR_DPnavspiralP1.Parameter.Labels.Index.dyn) + 1;
shots_per_volumn = profiles / ch_nr_p / n_nsa / n_dyn;
k_spa_P1_data = reshape(k_spa_P1_data,kx, ch_nr_p, n_nsa, shots_per_volumn, n_dyn);
[kx, n_ch,  n_nsa, shots, dynamics] = size(k_spa_P1_data)
k_spa_P1_data_rm_phase_offset = k_spa_P1_data;


MR_DPnavspiralP2 = MRecon(fn_P2);
MR_DPnavspiralP2.Parameter.Parameter2Read.typ = 1;
MR_DPnavspiralP2.Parameter.Parameter2Read.mix = 0;
MR_DPnavspiralP2.ReadData;
MR_DPnavspiralP2.RandomPhaseCorrection;
k_spa_P2_data = double(MR_DPnavspiralP2.Data);
[kx, profiles] = size(k_spa_P2_data);
ch_nr_p = length(MR_DPnavspiralP2.Parameter.Labels.CoilNrs);
n_nsa = max(MR_DPnavspiralP2.Parameter.Labels.Index.aver) + 1;
n_dyn = max(MR_DPnavspiralP2.Parameter.Labels.Index.dyn) + 1;
shots_per_volumn = profiles / ch_nr_p / n_nsa / n_dyn;
k_spa_P2_data = reshape(k_spa_P2_data,kx, ch_nr_p, n_nsa, shots_per_volumn, n_dyn);
[kx, n_ch,  n_nsa, shots, dynamics] = size(k_spa_P2_data)
k_spa_P2_data_rm_phase_offset = k_spa_P2_data;


for diffusion_nr = 1:dynamics
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

MR_DPnavspiralS1 = MRecon(fn_S1);
MR_DPnavspiralS1.Parameter.Parameter2Read.typ = 1;
MR_DPnavspiralS1.Parameter.Parameter2Read.mix = 0;
MR_DPnavspiralS1.ReadData;
MR_DPnavspiralS1.RandomPhaseCorrection;
k_spa_S1_data = double(MR_DPnavspiralS1.Data);
[kx, profiles] = size(k_spa_S1_data);
ch_nr_s = length(MR_DPnavspiralS1.Parameter.Labels.CoilNrs);
n_nsa = max(MR_DPnavspiralS1.Parameter.Labels.Index.aver) + 1;
n_dyn = max(MR_DPnavspiralS1.Parameter.Labels.Index.dyn) + 1;
shots_per_volumn = profiles / ch_nr_s / n_nsa / n_dyn;
k_spa_S1_data = reshape(k_spa_S1_data,kx, ch_nr_s, n_nsa, shots_per_volumn, n_dyn);
[kx, n_ch,  n_nsa, shots, dynamics] = size(k_spa_S1_data)
k_spa_S1_data_rm_phase_offset = k_spa_S1_data;


MR_DPnavspiralS2 = MRecon(fn_S2);
MR_DPnavspiralS2.Parameter.Parameter2Read.typ = 1;
MR_DPnavspiralS2.Parameter.Parameter2Read.mix = 0;
MR_DPnavspiralS2.ReadData;
MR_DPnavspiralS2.RandomPhaseCorrection;
k_spa_S2_data = double(MR_DPnavspiralS2.Data);
[kx, profiles] = size(k_spa_S2_data);
ch_nr_s = length(MR_DPnavspiralS2.Parameter.Labels.CoilNrs);
n_nsa = max(MR_DPnavspiralS2.Parameter.Labels.Index.aver) + 1;
n_dyn = max(MR_DPnavspiralS2.Parameter.Labels.Index.dyn) + 1;
shots_per_volumn = profiles / ch_nr_s / n_nsa / n_dyn;
k_spa_S2_data = reshape(k_spa_S2_data,kx, ch_nr_s, n_nsa, shots_per_volumn, n_dyn);
[kx, n_ch,  n_nsa, shots, dynamics] = size(k_spa_S2_data)
k_spa_S2_data_rm_phase_offset = k_spa_S2_data;

for diffusion_nr = 1:dynamics
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

sel_ch = [1:ch_nr_m];
% sel_ch =  [1:17 19:27]
for ch=1:length(sel_ch)
    im_M1_ksp_phase_unwrap(:,ch,:) = unwrap(squeeze(angle(im_M1_ksp(:,sel_ch(ch),:))));
    im_M2_ksp_phase_unwrap(:,ch,:) = unwrap(squeeze(angle(im_M2_ksp(:,sel_ch(ch),:))));
end

sel_ch = [1:ch_nr_p];
% sel_ch = [1 3:17 19:27]
for ch=1:length(sel_ch)
    im_P1_ksp_phase_unwrap(:,ch,:) = unwrap(squeeze(angle(im_P1_ksp(:,sel_ch(ch),:))));
    im_P2_ksp_phase_unwrap(:,ch,:) = unwrap(squeeze(angle(im_P2_ksp(:,sel_ch(ch),:))));
end

sel_ch = [1:ch_nr_s];
% sel_ch = [1 3:24 27]
for ch=1:length(sel_ch)
    im_S1_ksp_phase_unwrap(:,ch,:) = unwrap(squeeze(angle(im_S1_ksp(:,sel_ch(ch),:))));
    im_S2_ksp_phase_unwrap(:,ch,:) = unwrap(squeeze(angle(im_S2_ksp(:,sel_ch(ch),:))));
end

diffusion_nr = input(['diffusion setting nr (1-',num2str(dynamics),'): ']);
b0_diffusion_idx = 1;

all_wrapped = false;
while(~all_wrapped)
    
    mm = im_M1_ksp_phase_unwrap(:,:,diffusion_nr)-im_M2_ksp_phase_unwrap(:,:,b0_diffusion_idx);
    figure(1001); set(gcf,'DefaultAxesFontSize',8);
    plot(mm(1,:)); set(gca,'xtick',[1:2:40]); grid on;
    pp =im_P1_ksp_phase_unwrap(:,:,diffusion_nr)-im_P2_ksp_phase_unwrap(:,:,b0_diffusion_idx);
    figure(1002); set(gcf,'DefaultAxesFontSize',8);
    plot(pp(1,:)); set(gca,'xtick',[1:2:40]); grid on;
    ss =im_S1_ksp_phase_unwrap(:,:,diffusion_nr)-im_S2_ksp_phase_unwrap(:,:,b0_diffusion_idx);
    figure(1003);  set(gcf,'DefaultAxesFontSize',8);
    plot(ss(1,:)); set(gca,'xtick',[1:2:40]); grid on;
    
    all_wrapped = input('all phase properly warpped?(false/true): ');
    
    if(~all_wrapped)
        wrap_idx = input('wrap index for kx (phase-2*pi): ');
        im_M1_ksp_phase_unwrap(:,wrap_idx,diffusion_nr) = im_M1_ksp_phase_unwrap(:,wrap_idx,diffusion_nr)-2* pi;
        wrap_idx_p =input('wrap index for ky (phase-2*pi): ');
        im_P1_ksp_phase_unwrap(:,wrap_idx_p,diffusion_nr) = im_P1_ksp_phase_unwrap(:,wrap_idx_p,diffusion_nr) - 2* pi;
        wrap_idx_z = input('wrap index for kz (phase-2*pi): ');
        im_S1_ksp_phase_unwrap(:,wrap_idx_z,diffusion_nr) = im_S1_ksp_phase_unwrap(:,wrap_idx_z,diffusion_nr) - 2* pi;
    end
    
end

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

for diffusion_nr = 1:dynamics
    phi_M(:,diffusion_nr) = squeeze(mean(im_M1_ksp_phase_unwrap(:,:,diffusion_nr)-im_M2_ksp_phase_unwrap(:,:,b0_diffusion_idx),2));
    phi_P(:,diffusion_nr) = squeeze(mean(im_P1_ksp_phase_unwrap(:,:,diffusion_nr)-im_P2_ksp_phase_unwrap(:,:,b0_diffusion_idx),2));
    phi_S(:,diffusion_nr) = squeeze(mean(im_S1_ksp_phase_unwrap(:,:,diffusion_nr)-im_S2_ksp_phase_unwrap(:,:,b0_diffusion_idx),2));
end

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

% weight = sqrt(diff(squeeze(trj_meas_kx(:,diffusion_nr))).^2 + diff(squeeze(trj_meas_ky(:,diffusion_nr))).^2 + diff(squeeze(trj_meas_kz(:,diffusion_nr))).^2);

%%

if(exist(trj_save_fn)>0)
    save(trj_save_fn, 'trj_meas_kx','trj_meas_ky','trj_meas_kz','-append');
else
    save(trj_save_fn, 'trj_meas_kx','trj_meas_ky','trj_meas_kz');
end

end