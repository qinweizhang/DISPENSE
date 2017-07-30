%% LOAD DATA
clear
fn = 'dp_16072017_1621562_9_2_wip_dpsti_sad_2dV4.raw';
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

%% load goal object and parameters and waveforms
nav_kz_steps = MR_TSEDPnav_data_recon1.Parameter.GetParameter('EX_T2PREP_DPnav_zRes');

GR_nav_mc_0 =  MR_TSEDPnav_data_recon1.Parameter.GetObject('GR`t2prep_dp_Nav_mc_0_'); %pre phase wrap
GR_nav_0 =  MR_TSEDPnav_data_recon1.Parameter.GetObject('GR`t2prep_dp_Nav_0_'); %in m(x)
GR_nav_1 =  MR_TSEDPnav_data_recon1.Parameter.GetObject('GR`t2prep_dp_Nav_1_'); %in p(y)
GR_nav_2 =  MR_TSEDPnav_data_recon1.Parameter.GetObject('GR`t2prep_dp_Nav_2_'); %in s(z)
GR_nav_mc_2 =  MR_TSEDPnav_data_recon1.Parameter.GetObject('GR`t2prep_dp_Nav_mc_2_'); %pre phase wrap
AQ_dpnav = MR_TSEDPnav_data_recon1.Parameter.GetObject('AQ`DPNav'); 
% load waveforms
GR_spiral = load('spiral_GR.dat');
x_shape = GR_spiral(:,1) .* GR_nav_0.GetValue('wf_str');
y_shape = GR_spiral(:,2) .* GR_nav_1.GetValue('wf_str');
x_shape_invert = GR_spiral(:,3) .* GR_nav_0.GetValue('wf_str');
y_shape_invert = GR_spiral(:,4) .* GR_nav_1.GetValue('wf_str');

GR_dwell = 0.0064;   %in ms
calc_dwell = 0.0001; %in ms

kz_steps = nav_kz_steps.values;
if(mod(kz_steps, 2)==0) %even
    gr_x = repmat(cat(1, x_shape, x_shape_invert) , kz_steps/2, 1);
    gr_y = repmat(cat(1, y_shape, y_shape_invert) , kz_steps/2, 1);
else %odd
    gr_x = cat(1, repmat(cat(1, x_shape, x_shape_invert) , (kz_steps-1)/2, 1), x_shape);
    gr_y = cat(1, repmat(cat(1, y_shape, y_shape_invert) , (kz_steps-1)/2, 1), y_shape);
end

gr_x_fine = col(repmat(gr_x, 1, 64)');
gr_y_fine = col(repmat(gr_y, 1, 64)');

mc_0_str =  GR_nav_mc_0.GetValue('str');
mc_0_slope =  GR_nav_mc_0.GetValue('slope');
mc_0_lenc =  GR_nav_mc_0.GetValue('lenc');
mc_0 = cat(2, 0, linspace(0, mc_0_str, mc_0_slope./calc_dwell),linspace(mc_0_str, mc_0_str, mc_0_lenc./calc_dwell),linspace(mc_0_str, 0, mc_0_slope./calc_dwell));

mc_2_str =  GR_nav_mc_2.GetValue('str');
mc_2_slope =  GR_nav_mc_2.GetValue('slope');
mc_2_lenc =  GR_nav_mc_2.GetValue('lenc');
mc_2 = cat(2, 0, linspace(0, mc_2_str, mc_2_slope./calc_dwell),linspace(mc_2_str, mc_2_str, mc_2_lenc./calc_dwell),linspace(mc_2_str, 0, mc_2_slope./calc_dwell));

nav_2_str =  GR_nav_2.GetValue('str');
nav_2_slope =  GR_nav_2.GetValue('slope');
nav_2_lenc =  GR_nav_2.GetValue('lenc');
nav_2 = cat(2, 0, linspace(0, nav_2_str, nav_2_slope./calc_dwell),linspace(nav_2_str, nav_2_str, nav_2_lenc./calc_dwell),linspace(nav_2_str, 0, nav_2_slope./calc_dwell));


mc_0_start_time = GR_nav_mc_0.GetValue('time') - GR_nav_mc_0.GetValue('ref'); %GR x
nav_1_start_time = GR_nav_1.GetValue('time') - GR_nav_1.GetValue('ref'); %GR x
mc_2_start_time = GR_nav_mc_2.GetValue('time') - GR_nav_mc_2.GetValue('ref'); %GR z
AQ_start_time = AQ_dpnav.GetValue('time') - AQ_dpnav.GetValue('ref'); %AQ

start_time = min([mc_0_start_time, nav_1_start_time, mc_2_start_time, AQ_start_time]);

gr_x_fine_all = cat(1, mc_0',  gr_x_fine);
gr_y_fine_all = gr_y_fine;

nav_2_start_time = GR_nav_2.GetValue('time') - GR_nav_2.GetValue('ref');
mc_2_end_time = GR_nav_mc_2.GetValue('time') - GR_nav_mc_2.GetValue('ref') + GR_nav_mc_2.GetValue('dur');
nav_2_interval = GR_nav_2.GetValue('interval');
nav_2_repeat = GR_nav_2.GetValue('repetitions');
gr_z_fine_all = cat(1, mc_2', zeros(round((nav_2_start_time - mc_2_end_time)./calc_dwell), 1 ), repmat(cat(1, nav_2', zeros(round(nav_2_interval./calc_dwell), 1)), nav_2_repeat ,1));


gr_x_fine_all_time_match = cat(1, zeros(round((mc_0_start_time - start_time)./calc_dwell) , 1), gr_x_fine_all );
gr_y_fine_all_time_match = cat(1, zeros(round((nav_1_start_time - start_time)./calc_dwell) , 1), gr_y_fine_all );
gr_z_fine_all_time_match = cat(1, zeros(round((mc_2_start_time - start_time)./calc_dwell) , 1), gr_z_fine_all );

%% trajectory calculation for all aq samples
traj_kx_fine = zeros(1, max([length(gr_x_fine_all_time_match), length(gr_y_fine_all_time_match), length(gr_z_fine_all_time_match)]));
traj_ky_fine = traj_kx_fine;
traj_kz_fine = traj_kx_fine;

gamma = 42.576e3; %Hz/mT
% k = sum(gamma*G*t), gamma in Hz/mT, G in mT/m, t in S; k in 1/m
for idx = 2:length(traj_kx_fine)
    if (idx > length(gr_x_fine_all_time_match))
        gr_x_temp = 0;
    else
        gr_x_temp = gr_x_fine_all_time_match(idx);        
    end
    
    if (idx > length(gr_y_fine_all_time_match))
        gr_y_temp = 0;
    else
        gr_y_temp = gr_y_fine_all_time_match(idx);    
    end
    
    if (idx > length(gr_z_fine_all_time_match))
        gr_z_temp = 0;
    else
        gr_z_temp = gr_z_fine_all_time_match(idx);    
    end
    
    traj_kx_fine(idx) = traj_kx_fine(idx - 1) + (gamma * gr_x_temp * (calc_dwell/1e3));
    traj_ky_fine(idx) = traj_ky_fine(idx - 1) + (gamma * gr_y_temp * (calc_dwell/1e3));
    traj_kz_fine(idx) = traj_kz_fine(idx - 1) + (gamma * gr_z_temp * (calc_dwell/1e3));
end

aq_start_time = AQ_dpnav.GetValue('time') - AQ_dpnav.GetValue('ref');
aq_samples = AQ_dpnav.GetValue('samples');
aq_interval = AQ_dpnav.GetValue('interval');

aq_samples_idx = round(([0:aq_samples-1] .* aq_interval + (aq_start_time - start_time)) / calc_dwell); 

trj_meas_kx = traj_kx_fine(aq_samples_idx)' / 1000;  %in 1/mm
trj_meas_ky = traj_ky_fine(aq_samples_idx)' / 1000;
trj_meas_kz = traj_kz_fine(aq_samples_idx)' / 1000;

figure(100); plot3(trj_meas_kx, trj_meas_ky, trj_meas_kz);
xlabel('kx m^-^1'); ylabel('ky m^-^1'); zlabel('kz m^-^1');

save traj_ideal.mat trj_meas_kx trj_meas_ky trj_meas_kz


%% CALCULATE TRAJECTORY: code from reconFrame

if isempty( GridderPars.GridOvsFactor )
    grid_ovs       = kx_ovs / 1;
else
    grid_ovs = [];
end

%---------------------------------------------------------------------
% Spiral fine-tuning (scanner specific)
channel_delay  = 0;
phase_offset   = 0;
lambda         = 3;

%---------------------------------------------------------------------
% Definition of k-space trajectory is found in mpuspiral__g.c in
% Philips pulse programming environment
acq_samples    = round(length(kx_range(1):kx_range(2) ) / kx_ovs);
no_turns       = acq_samples /(2*no_interleaves);                       % ok
lambda         = min(no_turns,lambda);                                  % /Ta
phi_top        = pi*acq_samples*sqrt(1+lambda)/(no_interleaves*no_samples);
A              = no_interleaves/(2*pi);
k              = zeros(no_samples,no_interleaves, 3,'single');

%--------------------------------------------------
%-------------------
% Calculate spiral trajectory
for interleave=0:no_interleaves-1
    phi_0 = 2*pi*interleave/no_interleaves+phase_offset;
    for sample=fix(channel_delay):no_samples-1
        samp = channel_delay - fix(channel_delay) + sample;
        phi_t = phi_top*samp/sqrt(1+lambda*samp/no_samples);
        k(sample+1,interleave+1,1) = A*phi_t*cos(phi_t-phi_0);
        k(sample+1,interleave+1,2) = A*phi_t*sin(phi_t-phi_0);
    end
end
k = reshape(k, size(k,1), size(k,2), 1, size(k,3) );
