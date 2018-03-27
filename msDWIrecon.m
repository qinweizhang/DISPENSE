% msDWIrecon: reconstructed multishot phase incoherent DWI images
%
%
% INPUT
% k_spa:              cpx phase inhoerent k space data in [ky kz n_coil n_shots]
% k_spa (for LRT):    cpx phase inhoerent k space data in [ky kz n_coil n_shots 2]; 1st for nav data; 2nd for image data
% sensemaps:          cpx sense maps in [ky kz n_coil]
% phase_error_maps:   cpx phase maps for every shot in [ky kz 1 n_shots]
% pars:               reconstruction parameters structure
%
% OUTPUT
% image_corrected:    reconstructed image
%
% (c) Qinwei Zhang (q.zhang@amc.uva.nl) 2017 @AMC Amsterdam



function image_corrected = msDWIrecon(kspa, sense_map, phase_error, pars)

if(~iscell(kspa))
    [ky_dim, kz_dim, nc, ns, n_type] = size(kspa);
    
    if (nc==1)
        assert(length(size(sense_map)) == 2 && sum([ky_dim kz_dim]==size(sense_map)) ==  2); %sense_map size check
    else
        assert(length(size(sense_map)) == 3 && sum([ky_dim kz_dim nc]==size(sense_map)) ==  3); %sense_map size check
    end
    
    
    if(ns > 1)
        if(strcmp(pars.method, 'CG_SENSE_I')) %phase_error only used for CG
            assert(length(size(phase_error)) == 4 && sum([ky_dim kz_dim 1 ns]==size(phase_error)) ==  4); %phase_error size check
        end
    else
        warning('This is a single-shot dataset!')
    end
else % kspa is a cell
    
    assert(strcmp(pars.method, 'LRT'));
    assert(pars.LRT.mix_trajectory == 1);
    disp('msDWIrecon is LRT with mixed trajectory: data checking ommitted!')
    [ns, n_type] = size(kspa);
end



%% algorithms
if strcmp(pars.method, 'CG_SENSE_I')
    %% IMAGE SPACE CG_SENSE
    warning('Perform CG_SENSE_I reconstruction...');
    assert(n_type==1);

    %preprocessing
    if(nc>1)
        sense_map = permute((normalize_sense_map(permute(sense_map,[4 1 2 3]))),[2 3 4 1]); %normalized in 4th dimension; so permute and permute back
    end
    phase_error = normalize_sense_map(phase_error); %miss use normalized_sense_map function to normalized phase_error map
    mask = abs(kspa) > 0;
    
    if((sum(abs(sense_map(:)))>0)&&(sum(abs(phase_error(:)))>0))
        %scale k space incase they are ununiformly sampled; As we are min(||kspace-EI||^2),  but now kspace and I have different energy,this is
        %not good.
        sample_times = sum(mask,4);
        kspa = bsxfun(@rdivide, kspa, sample_times);
        kspa(find(isnan(kspa)))=0; kspa(find(isinf(kspa)))=0;
        
        
        b = col(kspa);
        
        
        A=FPSoperator(sense_map, phase_error, [ky_dim kz_dim], nc, ns, mask);
        %         image_corrected = A'*b;  %direct inverse
        lamda =  pars.CG_SENSE_I.lamda;
        maxit =  pars.CG_SENSE_I.nit;
        tol =   pars.CG_SENSE_I.tol;
        
        image_corrected_unfiltered=regularizedReconstruction(A,b,@L2Norm,lamda,'maxit',maxit,'tol', tol);
        image_corrected = filt_perifiral_kspa(image_corrected_unfiltered, squeeze(sum(mask(:,:,1,:),4)), pars.trim_kspa_filter_mask_size);
    else
        image_corrected = zeros(ky_dim, kz_dim);
    end
elseif strcmp(pars.method, 'CG_SENSE_K')    
    %% KSPACE CG_SENSE TODO
    
    
    
elseif strcmp(pars.method, 'POCS_ICE') %update phase map iteratively
    %% POCS_ICE 
     warning('Perform POCS_ICE reconstruction...no external phase error information is used...');
    assert(n_type==1);

    trj = [];
    image_corrected = POCS_SENESE(kspa,trj, sense_map, pars);
    
    
elseif strcmp(pars.method, 'LRT') %recon in the LRT frame
    %% LRT 
    warning('Perform LRT reconstruction...no external phase error information is used...');
    assert(n_type==2);

    image_corrected_LRT=LRT_recon_msDWI(kspa,squeeze(sense_map),pars.LRT);
    figure(1000);
    subplot(221); montage(permute(squeeze(abs(image_corrected_LRT(:,:,1,:,1))),[1 2 4 3]),'displayrange',[]); title('LRT recon of nav. column (mag.)')
    subplot(222); montage(permute(squeeze(angle(image_corrected_LRT(:,:,1,:,1))),[1 2 4 3]),'displayrange',[-pi pi]);title('LRT recon of nav. column (phase.)')
    subplot(223); montage(permute(squeeze(abs(image_corrected_LRT(:,:,1,:,2))),[1 2 4 3]),'displayrange',[]); title('LRT recon of image column (mag.)')
    subplot(224); montage(permute(squeeze(angle(image_corrected_LRT(:,:,1,:,2))),[1 2 4 3]),'displayrange',[-pi pi]); title('LRT recon of image column (phase.)')
    
    warning('!!! LRT reconed non-b0 images are averaged to obtain final results !!!');
    image_corrected = image_corrected_LRT;
%     image_corrected = mean(squeeze(image_corrected_LRT(:,:,1,2:end,2)),3);
elseif strcmp(pars.method, 'LRT_implicit') %recon in the LRT_implicit frame
     %% LRT_implicit 
    warning('Perform LRT implicit reconstruction...');

    image_corrected_LRT=LRT_implicit_recon_msDWI(kspa,squeeze(sense_map),pars.LRT_implicit);
    
    figure(1000);
    subplot(221); montage(permute(squeeze(abs(image_corrected_LRT(:,:,1,:,1))),[1 2 4 3]),'displayrange',[]); title('LRT recon of nav. column (mag.)')
    subplot(222); montage(permute(squeeze(angle(image_corrected_LRT(:,:,1,:,1))),[1 2 4 3]),'displayrange',[-pi pi]);title('LRT recon of nav. column (phase.)')
    subplot(223); montage(permute(squeeze(abs(image_corrected_LRT(:,:,1,:,2))),[1 2 4 3]),'displayrange',[]); title('LRT recon of image column (mag.)')
    subplot(224); montage(permute(squeeze(angle(image_corrected_LRT(:,:,1,:,2))),[1 2 4 3]),'displayrange',[-pi pi]); title('LRT recon of image column (phase.)')
    
    warning('!!! LRT reconed non-b0 images are averaged to obtain final results !!!');
    image_corrected = image_corrected_LRT;
else
    error('recon method not recognized...')
end

image_corrected(find(isnan(image_corrected))) = 0; image_corrected(find(isinf(image_corrected))) = 0; 
end