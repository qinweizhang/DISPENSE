% msDWIrecon: reconstructed multishot phase incoherent DWI images
%
%
% INPUT
% k_spa:              phase inhoerent k space data in [kx ky kz n_coil n_shots]
% sensemaps:          sense maps in [kx ky kz n_coil]
% phase_error_maps:   phase maps for every shot in [kx ky kz 1 n_shots]
% pars:          reconstruction parameters structure
%
% OUTPUT
% image_corrected:    reconstructed image
%
% (c) Qinwei Zhang (q.zhang@amc.uva.nl) 2017 @AMC



function image_corrected = msDWIrecon(kspa, sense_map, phase_error, pars)

[ky_dim, kz_dim, nc, ns] = size(kspa);

assert(length(size(sense_map)) == 3 && sum([ky_dim kz_dim nc]==size(sense_map)) ==  3); %sense_map size check

if(ns > 1)
    assert(length(size(phase_error)) == 4 && sum([ky_dim kz_dim 1 ns]==size(phase_error)) ==  4); %phase_error size check
else
    warning('This is a single-shot dataset!')
end

%normalized maps
sense_map = squeeze(normalize_sense_map(permute(sense_map,[4 1 2 3]))); %normalized in 4th dimension
phase_error = normalize_sense_map(phase_error); %miss use normalized_sense_map function to normalized phase_error map

%% recon
if strcmp(pars.method, 'CG_CENSE')
    mask = abs(kspa) > 0;
    A=FPSoperator(sense_map, phase_error, [ky_dim kz_dim], nc, ns, mask);
    
    b = col(kspa);
    lamda =  pars.lamda;
    maxit =  pars.nit;
    
    image_corrected=regularizedReconstruction(A,b,@L2Norm,lamda,'maxit',maxit,'tol', 1e-10);
    
    
elseif strcmp(pars.method, 'POCS_CENSE')
    %TODO
elseif strcmp(pars.method, 'LRT')
    %TODO
else
    error('recon method not recognized...')
end


end