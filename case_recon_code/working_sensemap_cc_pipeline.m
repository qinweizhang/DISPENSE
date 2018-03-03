%works for
%cd('/home/qzhang/lood_storage/divi/Users/qzhang/Data/2018_01_23_SND_brain')
%dataset
%%
clear; clc; close all
cd('/home/qzhang/lood_storage/divi/Users/qzhang/Data/2018_02_28_SND_brain_jasper')

%% Set path
clear; close all; clc

data_fn = 'br_28022018_1634494_2_2_wip_sc2_3d_snd_brain_4bV4.raw';
sense_ref_fn = 'br_28022018_1634181_1000_5_wip_senserefscanV4.raw';
coil_survey_fn  = 'br_28022018_1632271_1000_2_wip_coilsurveyscanV4.raw';


%% TSE get VirtualCoilMartix
close all; clc
disp(' TSE data sorting and default recon...')

parameter2read.dyn = [];
parameter2read.cc_nr = 4; %0 for no cc
parameter2read.sense_recon = 0;

[ima_k_spa_data,TSE.ky_matched,TSE.kz_matched,TSE.shot_matched, TSE.ch_dim,ima_kspa_sorted, ima_default_recon, TSE_sense_map, TSE.kxrange, TSE.kyrange, TSE.kzrange, TSE.VirtualCoilMartix] = ...
    TSE_data_sortting(data_fn, sense_ref_fn, coil_survey_fn,parameter2read);

figure(610); immontage4D(permute(abs(ima_default_recon(:,:,:,:)),[1 2 4 3]), []);

TSE
assert(length(TSE.ky_matched)==size(ima_k_spa_data,2),'Profile number does not match with data size!')
disp('-finished- ');

%% get SENSE map and compress it

TSE.SENSE_kx =1;
TSE.SENSE_ky =2;
TSE.SENSE_kz =1;

% TSE.kyrange = [-62 -1]; 
TSE.kxrange = [-352 -1]; %consider now the ima_k_spa_data is oversampled in kx; kx oversmapled by 2 + 
TSE.kzrange = [-68, -1];

TSE.Ixrange = [ceil(TSE.kxrange(1).*TSE.SENSE_kx) -1];
TSE.Iyrange = [ceil(TSE.kyrange(1).*TSE.SENSE_ky) -1];
TSE.Izrange = [ceil(TSE.kzrange(1).*TSE.SENSE_kz) -1];

os = [1, 1, 1];
dim = [range(TSE.Ixrange), range(TSE.Iyrange), range(TSE.Izrange) ]+1;

clc_sens_par.cc = 1;
clc_sens_par.disp = 0;
[sense_map_temp, TSE.sense_Psi] = get_sense_map_external(sense_ref_fn, data_fn, coil_survey_fn, [dim(1)/2 dim(2) dim(3)], os, clc_sens_par);
%----compress sense map and sense_Psi
if(isfield(TSE, 'VirtualCoilMartix'))
    if(~isempty(TSE.VirtualCoilMartix))
        [sense_map_temp_cc, sense_Psi_cc] = compress_sense_map_Psi(TSE.VirtualCoilMartix, sense_map_temp,  TSE.sense_Psi);
    end
end
        

figure; montage(permute(squeeze(abs(sense_map_temp(:,:,32,:))),[1 2 4 3]),'displayrange',[]); title('sense map before cc')
figure; montage(permute(squeeze(abs(sense_map_temp_cc(:,:,32,:))),[1 2 4 3]),'displayrange',[]); title('sense map after cc')
