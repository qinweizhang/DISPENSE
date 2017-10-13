function  [sens_map, varargout] = get_sense_map_external(sense_ref, data_fn, coil_survey, varargin)

if(nargin == 4)
    recon_dim = varargin{1};
end
if(nargin == 5)
    recon_dim = varargin{1};
    os = varargin{2};
else
    os = [1 1 1];
end

s = MRecon(sense_ref);
survey = MRecon(coil_survey);
r = MRecon(data_fn);

ori_KxOversampling = r.Parameter.Encoding.KxOversampling;
ori_KyOversampling = r.Parameter.Encoding.KyOversampling;
ori_Offcentre = r.Parameter.Scan.Offcentre;

r.Parameter.Encoding.KxOversampling = [os(1); os(1)]; %no oversampling in ky actually [check this every time]
r.Parameter.Encoding.KyOversampling = [os(2); os(2)]; %no oversampling in ky actually [check this every time]


new_Offcentre = r.Parameter.Scan.Offcentre;

warning_msg = sprintf('Kerry SENSE maps calc. : manually changed Kx, KyOversampling %1.2f %1.2f>> %1.2f %1.2f', ori_KxOversampling(1), ori_KyOversampling(1), os(1), os(2));
warning(warning_msg);
warning_msg = sprintf('Kerry SENSE maps calc. : manually offcenter %1.2f %1.2f %1.2f>> %1.2f %1.2f %1.2f', ori_Offcentre(1), ori_Offcentre(2), ori_Offcentre(3), new_Offcentre(1), new_Offcentre(2), new_Offcentre(3));
warning(warning_msg);
str = sprintf('FOV  : %1.2f %1.2f %1.2f', r.Parameter.Scan.FOV(1), r.Parameter.Scan.FOV(2), r.Parameter.Scan.FOV(3));
disp(str);
str = sprintf('SENSE: %1.2f %1.2f %1.2f', r.Parameter.Scan.SENSEFactor(1), r.Parameter.Scan.SENSEFactor(2), r.Parameter.Scan.SENSEFactor(3));
disp(str);


MR_sense = MRsense(s, r, survey);
MR_sense.Mask = 1;
MR_sense.Smooth = 1;
MR_sense.Extrapolate = 1;
MR_sense.MatchTargetSize = 1;
% match the size. not sure if this is the correct way to do so.???
rs = 0;
if(0) %(recon_dim(3) == 1) %at least 3 slices? for some case this is not true?
    rs = 1;
    recon_dim(3) = 3;
end
if(exist('recon_dim'))
    MR_sense.OutputSizeReformated = recon_dim;
    MR_sense.OutputSizeSensitivity = recon_dim;
end

MR_sense.Perform;

sens_map_default_sorting = double(MR_sense.Sensitivity);
sens_Psi_default_sorting = MR_sense.Psi;

sens_map_no_sorting = sens_map_default_sorting .* 0;
sens_Psi_no_sorting = sens_Psi_default_sorting .* 0;

ch_order = r.Parameter.Labels.CoilNrs(:,1); %targeting sorting order. 
sorted_order = sort(ch_order);  %current sorting order

for c =1:length(ch_order)
    permute_order(c) =  find(sorted_order==ch_order(c));
end

sens_map_no_sorting = sens_map_default_sorting(:,:,:,permute_order);
sens_Psi_no_sorting = sens_Psi_default_sorting(permute_order, permute_order);

% for c =1:length(ch_order)
%     ch_nr = sorted_order(c);
%     permute_oder = find(ch_order==ch_nr);
%     sens_map_no_sorting(:,:,:,permute_oder) = sens_map_default_sorting(:,:,:,c);
% end

sens_map = sens_map_no_sorting;

ch_order'
str = sprintf('Kerry: Self sense calculation sort channel dimention based on the natural channel order in list file/chan_labels/Parameter.Labels.CoilNrs. \n\t\ti.e. The channel label value is not considered. \n\t\t Should do the same for the data!  \n\t\t Psi is also reordered accordingly');
warning(str);

if(nargout==2)
    varargout{1} = MR_sense.Psi;
end
% sens_map = sens_map_default_sorting;

if(rs == 1)
    sens_map = sens_map(:,:,2,:);
end

% display_bool = input('display sense map? (true/false): ');
display_bool = 0;
if(display_bool)
    figure(711);
    subplot(211);
    immontage4D(abs(sens_map)); title('sens maps (mag)'); xlabel('channels'); ylabel('slices');
    subplot(212);
    immontage4D(angle(sens_map),[-pi pi]); title('sens maps (phase)'); xlabel('channels'); ylabel('slices');
    
    figure(712);
    subplot(121);
    imagesc(abs(sens_Psi_default_sorting)); title('sens Psi (defualt)');
    subplot(122);
    imagesc(abs(sens_Psi_no_sorting)); title('sens Psi (modifiled order)');
end

end