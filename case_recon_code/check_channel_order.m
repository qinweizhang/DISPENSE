
%% recon b0 image without cc and sense unfolding
mrtest_raw = MRecon(data_fn);

mrtest = mrtest_raw.Copy;
mrtest.Parameter.Parameter2Read.typ = 1;
mrtest.Parameter.Parameter2Read.mix = 0;
mrtest.Parameter.Parameter2Read.dyn = 0;
mrtest.Parameter.Recon.ImmediateAveraging = 'No';
mrtest.ReadData;
mrtest.RandomPhaseCorrection;
% MR_DPstiRemoveOversampling;
mrtest.PDACorrection;
mrtest.DcOffsetCorrection;
mrtest.MeasPhaseCorrection;
mrtest.SortData;
kspa = double(mrtest.Data);
for ii = 1:size(kspa, 2)
    for jj = 1:size(kspa, 3)
        kspa(:,ii,jj,:) =  kspa(:,ii,jj,:) .* (-1)^(ii+jj);
    end
end
image = ifft3d(kspa);

image_pad = cat(2, zeros(size(kspa,1), 30, size(kspa,3), size(kspa,4)), image, zeros(size(kspa,1), 30, size(kspa,3), size(kspa,4)));

figure; montage(permute(squeeze(abs((image_pad(:,:,20,:)))),[1 2 4 3]),'displayrange',[]);

mrtest.Parameter.Encoding.XRange(1,:)
mrtest.Parameter.Encoding.YRange(1,:)
mrtest.Parameter.Encoding.ZRange(1,:)
mrtest.Parameter.Labels.CoilNrs(:,1)'
length(mrtest.Parameter.Labels.CoilNrs(:,1)')
%% get an ACM for cc
mrACM = MRecon(data_fn);
mrACM.Parameter.Parameter2Read.typ = 1;
mrACM.Parameter.Parameter2Read.mix = 0;
mrACM.Parameter.Parameter2Read.dyn = 0;
mrACM.Parameter.Recon.ImmediateAveraging = 'No';
mrACM.Parameter.Recon.ArrayCompression='Yes';
mrACM.Parameter.Recon.ACNrVirtualChannels = 4;
ACM = mrACM.Parameter.Recon.ACMatrix;
%% get default sense map for cc

SENSE_kx =1;
SENSE_ky =2;
SENSE_kz =1;
kyrange = [-62 -1]; 
kxrange = [-352 -1]; %consider now the ima_k_spa_data is oversampled in kx; kx oversmapled by 2 + 
kzrange = [-68, -1];
Ixrange = [ceil(kxrange(1).*SENSE_kx) -1];
Iyrange = [ceil(kyrange(1).*SENSE_ky) -1];
Izrange = [ceil(kzrange(1).*SENSE_kz) -1];


os = [1, 1, 1];
dim = [range(Ixrange), range(Iyrange), range(Izrange) ]+1;

clc_sens_par.cc = 1;
clc_sens_par.disp = 0;
[sense_map_default, sense_Psi_default] = get_sense_map_external(sense_ref_fn, data_fn, coil_survey_fn, [dim(1)/2 dim(2) dim(3)], os, clc_sens_par);



%% check if match and find he correct calibration order

sensemap_calib_order = [1:26 28:32];
sense_map_reorder = sense_map_default(:,:,:,sensemap_calib_order);
sense_Psi_reorder = sense_Psi_default(sensemap_calib_order, sensemap_calib_order);

s_id = 11;
figure(1101); montage(permute(squeeze(abs((image_pad(80:270,:,s_id,:)))),[1 2 4 3]),'displayrange',[0 500],'size',[6 6 ]); title('DISPENSE b0 direct recon (no SENSE unfold)')
figure(1102); montage(permute(squeeze(abs((sense_map_default(:,:,s_id,:)))),[1 2 4 3]),'displayrange',[0 2],'size',[6 6 ]); title('defualt sense maps')
figure(1103); montage(permute(squeeze(abs((sense_map_reorder(:,:,s_id,:)))),[1 2 4 3]),'displayrange',[0 2],'size',[6 6 ]); title('reordered sense maps')


[sense_map_default_cc, sense_Psi_default_cc] = compress_sense_map_Psi(ACM, sense_map_default,  sense_Psi_default);
[sense_map_reorder_cc, sense_Psi_reorder_cc] = compress_sense_map_Psi(ACM, sense_map_reorder,  sense_Psi_reorder);

figure(1104);
subplot(121); montage(permute(squeeze(abs((sense_map_default_cc(:,:,s_id,:)))),[1 2 4 3]),'displayrange',[]); title('defualt sense maps cc')
subplot(122); montage(permute(squeeze(abs((sense_map_reorder_cc(:,:,s_id,:)))),[1 2 4 3]),'displayrange',[]); title('reordered sense maps cc')

