function nav_k_spa_data = nav_kspa_data_read(fn)
%
%OUT
%
%nav_k_spa_data:    raw navigator k-space data

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



end