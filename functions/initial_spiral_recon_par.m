function recon_par = initial_spiral_recon_par


recon_par.ignore_kz = 0;
recon_par.acq_dim = [42 42 10];
recon_par.recon_dim  = [42 42 10];
recon_par.dyn_nr = 1;
recon_par.skip_point = 0 ;
recon_par.end_point = []; %or []: till the end;
recon_par.selected_point = [];  %overrule skip_point and end_point
recon_par.interations = 10;
recon_par.lamda = 0;
recon_par.recon_all_shot = 0;
recon_par.sense_map_recon =1;
recon_par.update_SENSE_map = 0;
recon_par.sense_calc_method = 'external'; %'ecalib' or 'external'
recon_par.sense_os = [1 1];  %oversampling in x and y: control sense FOV
recon_par.data_fn = []; %data_fn;
recon_par.sense_ref = []; %sense_ref_fn;
recon_par.coil_survey = []; %coil_survey_fn;


recon_par.channel_by_channel = 1;
recon_par.channel_by_channel = recon_par.channel_by_channel .* (1-recon_par.sense_map_recon );

recon_par.update_SENSE_map = 0;


recon_par.time_segmented_recon_for_B0_inhomo = 1; %time segmented recon to compensate B0 inhomongneity; !!!time consuming!!!
recon_par.time_segments.nr_segments = 10;  %more time segments more acurate; but recon times increases by a factor of time_segments
recon_par.time_segments.aq_interval = 0.01; %delta_t in ms


end