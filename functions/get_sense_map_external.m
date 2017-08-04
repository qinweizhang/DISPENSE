function  sens_map = get_sense_map_external(sense_ref, data_fn, coil_survey, recon_dim)

s = MRecon(sense_ref);
survey = MRecon(coil_survey);
r = MRecon(data_fn);

ori_KyOversampling = r.Parameter.Encoding.KyOversampling;
ori_Offcentre = r.Parameter.Scan.Offcentre;

r.Parameter.Encoding.KyOversampling = [1; 1]; %no oversampling in ky actually [check this every time]
r.Parameter.Scan.Offcentre = [0 0 0];  % do not apply offcenter for sens map calculation, as spiral image always reconed at the isocenter

new_KyOversampling = r.Parameter.Encoding.KyOversampling;
new_Offcentre = r.Parameter.Scan.Offcentre;

warning_msg = sprintf('Kerry SENSE maps calc. : manually changed KyOversampling %f %f>> %f %f', ori_KyOversampling(1), ori_KyOversampling(2), new_KyOversampling(1), new_KyOversampling(2));
warning(warning_msg);
warning_msg = sprintf('Kerry SENSE maps calc. : manually offcenter %f %f %f>> %f %f %f', ori_Offcentre(1), ori_Offcentre(2), ori_Offcentre(3), new_Offcentre(1), new_Offcentre(2), new_Offcentre(3));
warning(warning_msg);

MR_sense = MRsense(s, r, survey);
MR_sense.Mask = 1;
MR_sense.MatchTargetSize = 1;
MR_sense.Perform;



%now sense map matches the size of TSE image
sens_map_tse = MR_sense.Sensitivity;

sens_map_tse_kspa = bart('fft 7', sens_map_tse);

bart_command_1 = sprintf('resize -c 0 %d 1 %d 2 %d',recon_dim(1),recon_dim(2),recon_dim(3) )
sens_map_spial_kspa = bart(bart_command_1, sens_map_tse_kspa);
sens_map = bart('fft -i 7', sens_map_spial_kspa);

figure(711); montage(permute(squeeze(abs(sens_map(:,:,:,1))),[1 2 4 3]),'displayrange',[]); title('sens maps all slice, 1 channel')

end