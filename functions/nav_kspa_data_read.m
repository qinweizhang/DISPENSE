function nav_k_spa_data = nav_kspa_data_read(fn)
%
%OUT
%
%nav_k_spa_data:    raw navigator k-space data

MR_TSEDPnav_data = MRecon(fn);

MR_TSEDPnav_data.Parameter.Parameter2Read.typ = 1;
MR_TSEDPnav_data.Parameter.Parameter2Read.mix = 1;  %for DPnav Spirals

MR_TSEDPnav_data.Parameter.Recon.ImmediateAveraging = 'No';
MR_TSEDPnav_data.ReadData;
MR_TSEDPnav_data.RandomPhaseCorrection;
MR_TSEDPnav_data.DcOffsetCorrection;
MR_TSEDPnav_data.MeasPhaseCorrection;

nav_k_spa_data = double(MR_TSEDPnav_data.Data);
[kx, profiles] = size(nav_k_spa_data);

ch_nr = length(MR_TSEDPnav_data.Parameter.Labels.CoilNrs);
n_nsa = max(MR_TSEDPnav_data.Parameter.Labels.Index.aver) + 1;
n_dyn = max(MR_TSEDPnav_data.Parameter.Labels.Index.dyn) + 1;
shots_per_volumn = profiles / ch_nr  / n_dyn;

nav_k_spa_data = reshape(nav_k_spa_data,kx, ch_nr, shots_per_volumn, n_dyn); %for normal Spirals

nav_k_spa_data = nav_k_spa_data(kx/2+1:end,:,:,:,:); %for DPnav when no pi jump happens

[kx, n_ch, shots, diffusion_setting] = size(nav_k_spa_data);
disp('======================Kspace Size==========================')
disp(['Kx length       : ', num2str(kx)]);
disp(['Channels        : ', num2str(n_ch)]);
disp(['Shots per volume: ', num2str(shots)]);
disp(['Voumes          : ', num2str(diffusion_setting)]);

%% display relevant spiral parameters
try
    disp('======================SPIRAL Pars.=========================')
    GOAL_value = MR_TSEDPnav_data.Parameter.GetValue('EX_T2PREP_DPnav_xyRes');
    disp(['DPnav_xyRes     : ', num2str(GOAL_value)]);
    GOAL_value = MR_TSEDPnav_data.Parameter.GetValue('EX_T2PREP_DPnav_zRes');
    disp(['DPnav_zRes      : ', num2str(GOAL_value)]);
    GOAL_value = MR_TSEDPnav_data.Parameter.GetValue('EX_T2PREP_DPnav_Rxy');
    disp(['DPnav_Rxy       : ', num2str(GOAL_value)]);
    GOAL_value = MR_TSEDPnav_data.Parameter.GetValue('EX_T2PREP_DPnav_Rz');
    disp(['DPnav_Rz        : ', num2str(GOAL_value)]);
    GOAL_Par = MR_TSEDPnav_data.Parameter.GetObject('AQ`DPnav');
    disp(['AQ samples      : ', num2str(GOAL_Par.GetValue('samples'))]);
    disp('====================Geom. Pars.============================')
    GOAL_value_1 = MR_TSEDPnav_data.Parameter.GetValue('EX_GEO_voi_orientation');
    GOAL_value_2 = MR_TSEDPnav_data.Parameter.GetValue('EX_GEO_cur_stack_fat_shift_dir');
    disp(['Scan Orientation: ', GOAL_value_1, ';  Phase Encoding direction: ', GOAL_value_2]);
    GOAL_value_1 = MR_TSEDPnav_data.Parameter.GetValue('EX_GEO_fov');
    GOAL_value_2 = MR_TSEDPnav_data.Parameter.GetValue('EX_GEO_fov_p');
    GOAL_value_3 = MR_TSEDPnav_data.Parameter.GetValue('EX_GEO_cur_stack_fov_s');
    GOAL_value_4 = MR_TSEDPnav_data.Parameter.GetValue('EX_GEO_slice_oversample_factor');
    disp(['FOV             M: ', num2str(GOAL_value_1),  '  P:', num2str(GOAL_value_2),  '   S*os: ',num2str(GOAL_value_3),' * ',num2str(GOAL_value_4)]);  
    GOAL_value_1 = MR_TSEDPnav_data.Parameter.GetValue('EX_GEO_voi_ap_offcentre');
    GOAL_value_2 = MR_TSEDPnav_data.Parameter.GetValue('EX_GEO_voi_fh_offcentre');
    GOAL_value_3 = MR_TSEDPnav_data.Parameter.GetValue('EX_GEO_voi_lr_offcentre');
    disp(['Offset          AP: ', num2str(GOAL_value_1),'   FH: ',num2str(GOAL_value_2),'   LR: ',num2str(GOAL_value_3)]);
    GOAL_value_1 = MR_TSEDPnav_data.Parameter.GetValue('EX_GEO_voi_ap_angulation');
    GOAL_value_2 = MR_TSEDPnav_data.Parameter.GetValue('EX_GEO_voi_fh_angulation');
    GOAL_value_3 = MR_TSEDPnav_data.Parameter.GetValue('EX_GEO_voi_lr_angulation');
    disp(['Angulation      AP: ', num2str(GOAL_value_1),'   FH: ',num2str(GOAL_value_2),'   LR: ',num2str(GOAL_value_3)]);
    disp('===========================================================')
    
catch
    warning('Kerry:  Parameter display not succeed!');
    
end
end