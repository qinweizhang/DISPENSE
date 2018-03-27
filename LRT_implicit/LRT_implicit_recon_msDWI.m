function I_recon=LRT_implicit_recon_msDWI(kspa,sense_map,params)
% Implicit Low-Rank Tensor reconstruction 

% recon is run as: I_recon=LRT_implicit_recon_msDWI(kspace,sens,params)
% params is a struct- should be made with params_init

% inputs:
% kspace: undersampled kspace with dimensions (x,y,coils,shots)
% sens  : sens maps - size (x,y,coils) 
% params: params is a struct- should be made with params_init

% based on Mani et al. Multi-shot Sensitivity-Encoded Diffusion Data Recovery Using Structured Low-Rank Matrix Completion (MUSSELS)

% 2018 Q Zhang - AMC Amsterdam


% model: m = argmin{ l2_norm(A(m)-y)        +         lamda * nuclearNorm(D)         +             beta * l21_norm(Phi(I))  }


%============== dimension check=================%
[nx, ny, nc, ns] = size(kspa);
[nx_sens, ny_sens, nc_sens] = size(sense_map);
[nx_pe, ny_pe, ~, ns_pe] = size(phase_error);

assert((nx == nx_sens)&(ny == ny_sens)&(ny == ny_pe)&(nx == nx_pe), 'image dimensions shoud be consistent!')
assert((nc == nc_sens), 'sense channel shoud be consistent!')
assert((ns == ns_pe), 'shot number shoud be consistent!')


%==============preprocessing=====================%
sampling_mask = squeeze(abs(kspa(:,:,1,:))>0);
sense_map_norm = squeeze(normalize_sense_map(permute(sense_map, [1 2 4 3])));



%% algorithm 

A  = A_opt(sense_map_norm, [nx, ny], nc, ns, sampling_mask);

%use of A:
%m = A' * kspa(:);  m in size of [nx, ny, ns], kspa is col [nx*ny*nc*ns, 1]

% initialize
m_0 = A' * kspa(:);
beta = 0;
gamma = 0;

% ADMM below



