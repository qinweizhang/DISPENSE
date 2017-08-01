function ima_kspa_data(fn, save_fn)


MR_TSEDPima_data = MRecon(fn);

MR_TSEDPima_data_recon1 = MR_TSEDPima_data.Copy;

MR_TSEDPima_data_recon1.Parameter.Parameter2Read.typ = 1;
MR_TSEDPima_data_recon1.Parameter.Parameter2Read.mix = 0;  %for DPnav Spirals

MR_TSEDPima_data_recon1.Parameter.Recon.ImmediateAveraging = 'No';
MR_TSEDPima_data_recon1.ReadData;
MR_TSEDPima_data_recon1.RandomPhaseCorrection;
MR_TSEDPima_data_recon1.DcOffsetCorrection;
MR_TSEDPima_data_recon1.MeasPhaseCorrection;

ima_k_spa_data = double(MR_TSEDPima_data_recon1.Data);
[kx, profiles] = size(ima_k_spa_data);

ch_nr = length(MR_TSEDPima_data_recon1.Parameter.Labels.CoilNrs);
n_nsa = max(MR_TSEDPima_data_recon1.Parameter.Labels.Index.aver) + 1;
n_dyn = max(MR_TSEDPima_data_recon1.Parameter.Labels.Index.dyn) + 1;
shots_per_volumn = profiles / ch_nr  / n_dyn;

ima_k_spa_data = reshape(ima_k_spa_data,kx, ch_nr, shots_per_volumn, n_dyn); %for normal Spirals

[kx, n_ch, shots, diffusion_setting] = size(ima_k_spa_data)



%% DEFAULT RECON
MR_TSEDPima_data_recon1.SortData;
MR_TSEDPima_data_recon1.GridData;
MR_TSEDPima_data_recon1.RingingFilter;
MR_TSEDPima_data_recon1.ZeroFill;
MR_TSEDPima_data_recon1.K2I;

MR_TSEDPima_data_recon1.GridderNormalization;
MR_TSEDPima_data_recon1.SENSEUnfold;
% MR_TSEDPima_data_recon1.PartialFourier; %Cause OUT OF MEMORY
MR_TSEDPima_data_recon1.ConcomitantFieldCorrection;
MR_TSEDPima_data_recon1.DivideFlowSegments;
% MR_TSEDPima_data_recon1.Parameter.Recon.CoilCombination = 'pc'; MR_TSEDPima_data_recon1.CombineCoils;
% MR_TSE_recon1.GeometryCorrection; %don't do this here
% MR_TSEDPima_data_recon1.RemoveOversampling; %weird SENSE factor....
MR_TSEDPima_data_recon1.ZeroFill;
MR_TSEDPima_data_recon1.FlowPhaseCorrection;
MR_TSEDPima_data_recon1.RotateImage;

im_data = squeeze(MR_TSEDPima_data_recon1.Data);

MR_TSEDPima_data_recon1.Parameter.Recon.CoilCombination = 'pc'; MR_TSEDPima_data_recon1.CombineCoils;
MR_TSEDPima_data_recon1.ShowData;

im_data_cc = double(MR_TSEDPima_data_recon1.Data);

%% SAVE DATA

if(exist(save_fn)>0)
    save(save_fn, 'ima_k_spa_data','im_data_cc','-append');
else
    save(save_fn, 'ima_k_spa_data','im_data_cc');
end


end