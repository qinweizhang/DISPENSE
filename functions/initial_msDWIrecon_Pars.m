function pars = initial_msDWIrecon_Pars

%% General
pars.sense_map = 'external';
pars.data_fn = [];
pars.sense_ref = [];
pars.coil_survey = [];
pars.nav_phase_sm_kernel = 3;  %3 or 5, 1:no soomthing
pars.recon_x_locs = []; %80:270;
pars.enabled_ch = []; %1:TSE.ch_dim;
pars.b0_shots = []; %[] means first dynamic
pars.nonb0_shots = [];

pars.method='CG_SENSE_I';

%% CG_SENSEI

pars.CG_SENSE_I.lamda=0.1;
pars.CG_SENSE_I.nit=10;
pars.CG_SENSE_I.tol = 1e-10;

%% CG_SENSEK

pars.CG_SENSE_K.lamda=0.1;
pars.CG_SENSE_K.nit=10;
pars.CG_SENSE_K.tol = 1e-10;

%% POCS

pars.POCS.Wsize =[10 10];
pars.POCS.lamda = 0.1;
pars.POCS.nit  = 10;
pars.POCS.tol = 1e-10;
pars.POCS.nufft = false;

%% LRT

pars.LRT = LRT_params_init();
pars.LRT.NUFFT_nav_sense =[];
pars.LRT.NUFFT_nav_1ch = [];
pars.LRT.trj_length = [];

end



