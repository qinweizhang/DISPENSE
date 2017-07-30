MR_DPsti_TSE_b = MRecon('dp_24102016_1626539_7_2_wipdtidpstitsedualscansenseV4.raw');

MR_DPsti_TSE_b.Parameter.Parameter2Read.typ = 1;
MR_DPsti_TSE_b.Parameter.Parameter2Read.chan = 10;
MR_DPsti_TSE_b.ReadData;
MR_DPsti_TSE_b.RandomPhaseCorrection;
MR_DPsti_TSE_b.RemoveOversampling;
MR_DPsti_TSE_b.PDACorrection;
MR_DPsti_TSE_b.DcOffsetCorrection;
MR_DPsti_TSE_b.MeasPhaseCorrection;

% DP_Ks_data = double(MR_DPsti_TSE_b.Data);
% shot_useful = shot(36:end);
% shot_useful = shot_useful(1:15:end);
% j = sqrt(-1);
% for ix = 1:length(tran)
%     cor_temp = find(shot_useful == ix);
%     cor_idx = cor_temp(find(cor_temp<=length(shot_useful)/2));
%     DP_Ks_data(:,find(shot_useful == ix)) = DP_Ks_data(:,find(shot_useful == ix))...
%         .*exp(-j*tran(ix));
% end
% MR_DPsti_TSE_b.Data = single(DP_Ks_data);

MR_DPsti_TSE_b.SortData;
MR_DPsti_TSE_b.GridData;
MR_DPsti_TSE_b.RingingFilter;
MR_DPsti_TSE_b.ZeroFill;
MR_DPsti_TSE_b.K2IM;
MR_DPsti_TSE_b.EPIPhaseCorrection;
MR_DPsti_TSE_b.K2IP;
MR_DPsti_TSE_b.GridderNormalization;
MR_DPsti_TSE_b.SENSEUnfold;  %no sense here
MR_DPsti_TSE_b.PartialFourier;
MR_DPsti_TSE_b.ConcomitantFieldCorrection;
MR_DPsti_TSE_b.DivideFlowSegments;
MR_DPsti_TSE_b.Average;
MR_DPsti_TSE_b.GeometryCorrection;
MR_DPsti_TSE_b.RemoveOversampling;
MR_DPsti_TSE_b.ZeroFill;
MR_DPsti_TSE_b.FlowPhaseCorrection;
MR_DPsti_TSE_b.RotateImage;
MR_DPsti_TSE_b.ShowData;

%%
%ref
clear DP_uncorrected_d1

DP_uncorrected_d1 = bart('nufft -i -t', double((traj_matrix_ideal_rs)), ...
    double(DP_Ks_data_all_dy1_pm));
DP_uncorrected_d1 = fftshift(DP_uncorrected_d1);

figure; montage(permute(fftshift(fftshift(abs(DP_uncorrected_d1),3),1),[1 2 4 3]),'Displayrange',[])


%ref
clear DP_corrected_d1

DP_corrected_d1 = bart('nufft -i -t', double((traj_matrix_cor_rs)), ...
    double(DP_Ks_data_all_cor_pm));
DP_corrected_d1 = fftshift(DP_corrected_d1);

figure; montage(permute(fftshift(fftshift(abs(DP_corrected_d1),3),1),[1 2 4 3]),'Displayrange',[])


