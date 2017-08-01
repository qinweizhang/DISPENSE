function  sens_map = get_sense_map_external(sense_ref, data_fn, coil_survey, recon_dim)

MR_sense = MRsense(sense_ref, data_fn, coil_survey);
MR_sense.Mask = 1;
MR_sense.MatchTargetSize = 1;
MR_sense.Perform;


%now sense map matches the size of TSE image
sens_map_tse = MR_sense.Sensitivity;

sens_map_tse_kspa = bart('fft 7', sens_map_tse);

bart_command_1 = sprintf('resize -c 0 %d 1 %d 2 %d',recon_dim(1),recon_dim(2),recon_dim(3) )
sens_map_spial_kspa = bart(bart_command_1, sens_map_tse_kspa);
sens_map = bart('fft -i 7', sens_map_spial_kspa);


figure; montage(abs(sens_map(:,:,11,:)),'displayrange',[])

end