function [ima_k_spa_data,ky_matched,kz_matched,shot_matched, ch_dim,kspa_sorted, ima_data, TSE_sense_map, varargout] = TSE_data_sortting(raw_data_fn, sense_ref_fn, coil_survey_fn,varargin)

%INPUT
%
%raw_data_fn: .raw file name
%
%
%OUTPUT
%
%ima_k_spa_data:    unsorted TSE kspace data. To be used in the following correction steps
%ky_matched:        ky labels mached ima_k_spa_data
%kz_matched:        kz labels mached ima_k_spa_data
%shot_matched:      shot number labels mached ima_k_spa_data
%ch_dim:            channel numbers
%kspa_sorted:       sorted kspace
%ima_data:          uncorrectd image data.
%TSE_sense_map:     TSE sense maps for TSE recon
%varargout:         optional: kx_range ky_range kz_range

%% Presetting
MR_DPstiTSE = MRecon(raw_data_fn);


MR_DPstiTSE.Parameter.Parameter2Read.typ = 1;
MR_DPstiTSE.Parameter.Parameter2Read.mix = 0;
MR_DPstiTSE.Parameter.Recon.ImmediateAveraging = 'No';

if(nargin==4)
    parameter2read = varargin{1};
    if(~isempty(parameter2read.dyn))
        MR_DPstiTSE.Parameter.Parameter2Read.dyn = parameter2read.dyn;
    end
end

    
ch_dim = length(MR_DPstiTSE.Parameter.Labels.CoilNrs);

if(nargout == 12)
    % ask for add coil compression
    virtual_coil_nr = input(['Please input virtual coil nr for TSE DATA (0~',num2str(ch_dim),'; 0 for no coil comprsion):']);
    
    if(virtual_coil_nr>0)
        MR_DPstiTSE.Parameter.Recon.ArrayCompression='Yes';
        MR_DPstiTSE.Parameter.Recon.ACNrVirtualChannels=virtual_coil_nr;
        ch_dim = virtual_coil_nr;
    end
end


%% Read Data


MR_DPstiTSE.ReadData;
MR_DPstiTSE.RandomPhaseCorrection;
% MR_DPstiTSE.RemoveOversampling;
MR_DPstiTSE.PDACorrection;
MR_DPstiTSE.DcOffsetCorrection;
MR_DPstiTSE.MeasPhaseCorrection;

%% Extract\Create Labels

%--------------creat a label for shot number--------------&
dyn=MR_DPstiTSE.Parameter.Labels.Index.dyn;
typ=MR_DPstiTSE.Parameter.Labels.Index.typ;
mix=MR_DPstiTSE.Parameter.Labels.Index.mix;
nsa=MR_DPstiTSE.Parameter.Labels.Index.aver;
ch=MR_DPstiTSE.Parameter.Labels.Index.chan;
readdata_ix=find((typ==1)&(mix == 0)&ismember(ch,MR_DPstiTSE.Parameter.Parameter2Read.chan)&ismember(dyn,MR_DPstiTSE.Parameter.Parameter2Read.dyn));
type_1_ix = find(typ==1);
shot = zeros(length(typ),1);

rf=MR_DPstiTSE.Parameter.Labels.Index.rf;
first_echoes_ix = find(rf==0);
ch_num = MR_DPstiTSE.Parameter.Labels.CoilNrs; %total number of the channels
first_echoes_ty1_ix = first_echoes_ix(find(first_echoes_ix>=type_1_ix(1))); %first echoes for every shots
first_echoes_ty1_first_ch = first_echoes_ty1_ix(1:length(ch_num)*1:end); %first lables for every shots; *2 means nav data 14ch + ima data 14ch
%     first_echoes_ty1_first_ch = first_echoes_ty1_ix(1:length(ch_num)*1:end); %first lables for every shots; when no DPnav enabled

k_line_per_shot = first_echoes_ty1_first_ch(2) - first_echoes_ty1_first_ch(1);
for k=first_echoes_ty1_first_ch(1):length(shot)
    shot(k)=ceil((k-first_echoes_ty1_first_ch(1) + 1) / k_line_per_shot);
end
disp(['total shots: ', num2str(max(shot))]);

%-------------match labels ky kz shot to data
ky_all = MR_DPstiTSE.Parameter.Labels.Index.ky;
kz_all = MR_DPstiTSE.Parameter.Labels.Index.kz;

ky_matched = ky_all(readdata_ix);
kz_matched = kz_all(readdata_ix);
shot_matched = shot(readdata_ix);
nsa_matched = nsa(readdata_ix);

%----------comsider coil compression-------------
% ky_matched kz_matched shot_matched nsa_matched are match with full data  set (i.e. when all the k lines are exported)

if(virtual_coil_nr>0)
    ky_matched_rs = reshape(ky_matched, length(ch_num), length(ky_matched)/length(ch_num));
    kz_matched_rs = reshape(kz_matched, length(ch_num), length(kz_matched)/length(ch_num));
    shot_matched_rs = reshape(shot_matched, length(ch_num), length(shot_matched)/length(ch_num));
    nsa_matched_rs = reshape(nsa_matched, length(ch_num), length(nsa_matched)/length(ch_num));
    
    ky_matched_rs2 = repmat(ky_matched_rs(1,:), virtual_coil_nr, 1 );
    kz_matched_rs2 = repmat(kz_matched_rs(1,:), virtual_coil_nr, 1 );
    shot_matched_rs2 = repmat(shot_matched_rs(1,:), virtual_coil_nr, 1 );
    nsa_matched_rs2 = repmat(nsa_matched_rs(1,:), virtual_coil_nr, 1 );
    
    ky_matched = ky_matched_rs2(:);
    kz_matched = kz_matched_rs2(:);
    shot_matched = shot_matched_rs2(:);
    nsa_matched = nsa_matched_rs2(:);
end
%----------comsider coil compression END-------------





%============Display shot labels=============

%Sort Shot Label in the same way as soring data
MR_DPsti_shotLabel = MR_DPstiTSE.Copy;
Data_temp = MR_DPsti_shotLabel.Data;
[kx_dim, kx_lines] = size(Data_temp);
shot_label = shot(readdata_ix);
shot_label_data = repmat(shot_label', [kx_dim, 1]);
% size(shot_label_data) should have the same size as MR_DPstiTSE.Data
MR_DPsti_shotLabel.Data = single(shot_label_data);
MR_DPsti_shotLabel.SortData;
Ima_Sorted_shot_lable = squeeze(double(MR_DPsti_shotLabel.Data));
clear shot_label_temp; shot_label_temp(:,:,1,:) = squeeze(Ima_Sorted_shot_lable(32,:,:,2,:));
figure(605); montage(shot_label_temp,'displayrange',[]); colormap jet; title('shot labels');colorbar
figure(606); imagesc(squeeze(Ima_Sorted_shot_lable(32,:,:,1,1))'); colormap jet; colorbar; axis off; axis equal; xlabel('ky'); ylabel('kz')

%-------------------DONE---------------------------

% >>>>>>>exported kspace data from here
% ima_k_spa_data_temp = squeeze(double(MR_DPstiTSE.Data));
ima_k_spa_data = squeeze(double(MR_DPstiTSE.Data));
% for prof_nr = 1:length(ky_all)
%     ima_k_spa_data = ima_k_spa_data_temp(:,prof_nr).* (-1)^(ky_all(prof_nr) + kz_all(prof_nr));
% end



%% ==============Continue with Image Kspa data extracting
MR_DPstiTSE.SortData;
kspa_sorted = squeeze(double(MR_DPstiTSE.Data));
MR_DPstiTSE.GridData;
MR_DPstiTSE.RingingFilter;
MR_DPstiTSE.ZeroFill;
MR_DPstiTSE.K2IM;
MR_DPstiTSE.EPIPhaseCorrection;
MR_DPstiTSE.K2IP;
MR_DPstiTSE.GridderNormalization;

sense_recon = 0;
sense_recon = input('SENSE unfold? ');
if(sense_recon)
    if(virtual_coil_nr==0) %only for no coil compression
        % --------------Calculate SENSE object-------------
        x = size(MR_DPstiTSE.Data,1);
        y = size(MR_DPstiTSE.Data,2);
        z = size(MR_DPstiTSE.Data,3);
        MR_sense = MRsense(sense_ref_fn, raw_data_fn, coil_survey_fn);
        MR_sense.Mask = 1;
        MR_sense.MatchTargetSize = 1;
        MR_sense.OutputSizeReformated = [x y z]; warning('kerry temp. overwrite defualt parameter here!');
        MR_sense.OutputSizeSensitivity = [x y z]; warning('kerry temp. overwrite defualt parameter here!');
        MR_sense.Perform;
        % ----------------------end-----------------------
        MR_DPstiTSE.Parameter.Recon.Sensitivities = MR_sense;

        MR_DPstiTSE.SENSEUnfold; 
    else
        warning('Kerry: SENSE unfolder not available for compressed data!'); MR_DPstiTSE.CombineCoils;
    end
else
    warning('Kerry: SENSE unfolder disabled for this case!'); MR_DPstiTSE.CombineCoils;
end

MR_DPstiTSE.ConcomitantFieldCorrection;
MR_DPstiTSE.DivideFlowSegments;

% MR_DPstiTSE.ZeroFill;  
warning('Kerry: zeros filling unfolder disabled for this case!'); 
MR_DPstiTSE.FlowPhaseCorrection;
% MR_DPstiTSE.RotateImage;
MR_DPstiTSE.Average;
MR_DPstiTSE.ShowData;
ima_data = squeeze(double(MR_DPstiTSE.Data));
if(exist('MR_sense'))
    TSE_sense_map = MR_sense.Sensitivity;
    TSE_sense_map = TSE_sense_map./max(abs(TSE_sense_map(:)));
else
    TSE_sense_map =[];
end

if (nargout>=9)
    varargout{1} =  MR_DPstiTSE.Parameter.Encoding.XRange(1,:);
end
if(nargout>=10)
    varargout{2} = MR_DPstiTSE.Parameter.Encoding.YRange(1,:);
end
if(nargout>=11)
    if(isempty(MR_DPstiTSE.Parameter.Encoding.ZRange)) %2D scan maybe
        varargout{3} = [0 0];
    else
        varargout{3} = MR_DPstiTSE.Parameter.Encoding.ZRange(1,:);
    end
end
if(nargout>=12) %TSE_VirtualCoilMartix
    varargout{4} = MR_DPstiTSE.Parameter.Recon.ACMatrix;
end
        

end