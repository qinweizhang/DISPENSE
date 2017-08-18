function  sens_map = get_sense_map_external(sense_ref, data_fn, coil_survey, varargin)

if(nargin == 4)
    recon_dim = varargin{1};
end

s = MRecon(sense_ref);
survey = MRecon(coil_survey);
r = MRecon(data_fn);

ori_KyOversampling = r.Parameter.Encoding.KyOversampling;
ori_Offcentre = r.Parameter.Scan.Offcentre;

r.Parameter.Encoding.KyOversampling = [1; 1]; %no oversampling in ky actually [check this every time]
r.Parameter.Scan.Offcentre = [0 0 0];  % do not apply offcenter for sens map calculation, as spiral image always reconed at the isocenter


new_KyOversampling = r.Parameter.Encoding.KyOversampling;
new_Offcentre = r.Parameter.Scan.Offcentre;

warning_msg = sprintf('Kerry SENSE maps calc. : manually changed KyOversampling %1.2f %1.2f>> %1.2f %1.2f', ori_KyOversampling(1), ori_KyOversampling(2), new_KyOversampling(1), new_KyOversampling(2));
warning(warning_msg);
warning_msg = sprintf('Kerry SENSE maps calc. : manually offcenter %1.2f %1.2f %1.2f>> %1.2f %1.2f %1.2f', ori_Offcentre(1), ori_Offcentre(2), ori_Offcentre(3), new_Offcentre(1), new_Offcentre(2), new_Offcentre(3));
warning(warning_msg);


MR_sense = MRsense(s, r, survey);
MR_sense.Mask = 1;
MR_sense.Smooth = 1;
MR_sense.Extrapolate = 1;
MR_sense.MatchTargetSize = 1;
% match the size. not sure if this is the correct way to do so.???
if(exist('recon_dim'))
    MR_sense.OutputSizeReformated = recon_dim;
    MR_sense.OutputSizeSensitivity = recon_dim;
end

MR_sense.Perform;

sens_map = double(MR_sense.Sensitivity);

figure(711);
subplot(211);
immontage4D(abs(sens_map)); title('sens maps (mag)'); xlabel('channels'); ylabel('slices');
subplot(212);
immontage4D(angle(sens_map),[-pi pi]); title('sens maps (phase)'); xlabel('channels'); ylabel('slices');

end