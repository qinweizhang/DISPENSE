function pars = initial_msDWIrecon_Pars


pars.method='CG_SENSE_I';

pars.CG_SENSE_I.lamda=0.1;
pars.CG_SENSE_I.nit=10;
pars.CG_SENSE_I.tol = 1e-10;

pars.CG_SENSE_K.lamda=0.1;
pars.CG_SENSE_K.nit=10;
pars.CG_SENSE_K.tol = 1e-10;

pars.POCS.kernelsize =[20 20];
pars.POCS.lamda = 0.1;
pars.POCS.nit  = 10;
pars.POCS.tol = 1e-10;


end



