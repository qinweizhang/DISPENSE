function params= LRT_params_init()
% initializes relevant parameters for the LRT reconstruction 
% possible to change these before running recon

params.inspectLg=true;              % option to dynamically choose spatial rank based on first guess of C
params.increase_penalty_parameters=false; %option to increase alpha and beta inside the loop
params.scaleksp=true;              %option to scale kspace before recon - for consistent parameter use (TODO!)

params.Lg=1;                        %spatial rank
params.L3=3;                        %rank of first parameter dimension
params.L4=3;                        %rank of second parameter dimension

params.sparsity_transform='TV';     %'TV'/'wavelet'/'TVOP'
params.niter=20;                    %number of outer iterations in algo

%initialize parameters
params.alpha= 2;                    %penalty parameter >0
params.beta=  2;                    %penalty parameter >0
params.lambda=1e-2;                 %sparsity parameter
params.mu=1e1;                      %sparsity parameter

params.Imref=[];                    %possible reference (gold standard) image
params.x=20;                        %pixel to plot during recon loop
params.y=20;                        %pixel to plot during recon loop

params.C.tol=1e-10;
params.C.maxiter=50;
params.G.tol=1e-13;
params.G.maxiter=10;
params.G.precon=true;               %optional preconditioning of G

params.subspacedim1=1;               % dimension along which to take f.s. vals 
params.subspacedim2=1;               % dimension along which to take f.s. vals

params.nav_estimate_1=[];
params.nav_estimate_2=[];
params.eigenvals_1=[];
params.eigenvals_2=[];

params.autolambda=0   ; 
params.automu=0       ;             % automatically estimate mu on s.t. of first iter
params.normalize_sense=1;           %automatically normalizes sense maps 

params.mix_trajectory =0;           %in case of data consist of different trajectory. e.g. spiral + cartesian + radial...
