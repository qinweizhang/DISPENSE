%     Optimized DPsti_TSE_phase_error_cor function for large data recon. Avoid huge matrix allocation
%
%     Use msDWIrecon function to correct phase inherence during multishot DPsti-TSE
%
%     INPUT
%
%     ima_k_spa_data:             Unsorted TSE kspace data in [kx, nprofiles]
%     TSE:                        Structure contains labels (ky, kz, shot, channel, dyn) for every profiles in ima_k_spa_data
%     TSE_sense_map:              sense maps from sense reference scan
%     nav_data:                   navigator cpx images for every shot in [nav_x nav_y nav_z shots]
%
%     OUTPUT
%
%     image_corrected:            Corrected DPsti-TSE images
%
%     (c) Qinwei Zhang (q.zhang@amc.uva.nl) 2017 @AMC Amsterdam

% TODO make DPsti_TSE_phase_error_cor for POCS_ICE option
function image_corrected = DPsti_TSE_phase_error_cor_large_scale(ima_k_spa_data, TSE, TSE_sense_map, nav_data, pars)

%% Data check

max_shot = max(TSE.shot_matched);
assert(max_shot == size(nav_data, 4));
assert(TSE.ch_dim == size(TSE_sense_map, 4)||isempty(TSE_sense_map));

profiles_per_dyn = size(ima_k_spa_data,2)/TSE.dyn_dim;
assert(round(profiles_per_dyn) == profiles_per_dyn);
shots_per_dyn = max_shot/TSE.dyn_dim;
assert(round(shots_per_dyn) == shots_per_dyn);
assert(max(pars.enabled_ch)<=TSE.ch_dim&&min(pars.enabled_ch)>0);

TSE.kx_dim = TSE.kxrange(2) - TSE.kxrange(1) + 1;
assert(TSE.kx_dim >=  size(ima_k_spa_data,1));


TSE.ky_dim = TSE.kyrange(2) - TSE.kyrange(1) + 1;  %max_ky * 2 + 1;
TSE.kz_dim = TSE.kzrange(2) - TSE.kzrange(1) + 1;

TSE.Ix_dim = TSE.Ixrange(2) - TSE.Ixrange(1) + 1;
TSE.Iy_dim = TSE.Iyrange(2) - TSE.Iyrange(1) + 1;  %max_ky * 2 + 1;
TSE.Iz_dim = TSE.Izrange(2) - TSE.Izrange(1) + 1;

nonb0_shots_range = pars.nonb0_shots;

%% Preprocessing on sense data

disp('Preprocessing on sense data');
tic

% %sense mask || now it should be calculated outside and stored in TSE structure
if(isfield(TSE, 'sense_mask'))
    if(isempty(TSE.sense_mask)) %if empty calc again
        dim = [TSE.Ix_dim/2 TSE.Iy_dim TSE.Iz_dim];
        os = [1, 1, 1];
        sense_map_temp = get_sense_map_external(pars.sense_ref, pars.data_fn, pars.coil_survey, dim, os);
        rs_command = sprintf('resize -c 0 %d', TSE.kx_dim);
        sense_map_temp = bart(rs_command, sense_map_temp);
        
        TSE.sense_mask = abs(sense_map_temp(:,:,:,1 ))>0;
        clear sense_map_temp;
    end
else %if not exist, calc again
    dim = [TSE.Ix_dim/2 TSE.Iy_dim TSE.Iz_dim];
    os = [1, 1, 1];
    sense_map_temp = get_sense_map_external(pars.sense_ref, pars.data_fn, pars.coil_survey, dim, os);
    rs_command = sprintf('resize -c 0 %d', TSE.Ix_dim);
    sense_map_temp = bart(rs_command, sense_map_temp);
    
    TSE.sense_mask = abs(sense_map_temp(:,:,:,1 ))>0;
    clear sense_map_temp;
end

%sense maps
if(length(pars.enabled_ch)==1) %one channel
    warning('Kerry: This is one channel recon!')
    sense_map_3D = ones(TSE.Ix_dim,TSE.Iy_dim,TSE.Iz_dim);
    
    
else
    %estimate sense maps
    
    if(strcmp(pars.sense_map, 'ecalib'))
        
        ecalib_sense_map_3D = bart('ecalib -S -m1 -c0.2', kspa_b0);
        figure(3);
        displayslice = round(size(ecalib_sense_map_3D, 3)/2);
        subplot(211);montage(angle(ecalib_sense_map_3D(:,:,displayslice,:)),'displayrange',[-pi pi]); title('ecalib sense map (phase)')
        subplot(212);montage(abs(ecalib_sense_map_3D(:,:,displayslice,:)),'displayrange',[]); title('ecalib sense map (mag.)')
        ecalib_sense_map_3D = normalize_sense_map(ecalib_sense_map_3D);
        
        sense_map_3D = ecalib_sense_map_3D;
        
        
    elseif(strcmp(pars.sense_map, 'external'))
        if(isempty(TSE_sense_map))
            dim = [TSE.Ix_dim/2 TSE.Iy_dim TSE.Iz_dim];
            os = [1, 1, 1];
            [TSE_sense_map,TSE.sense_Psi]  = get_sense_map_external(pars.sense_ref, pars.data_fn, pars.coil_survey, dim, os);
            %compress sense map and sense_Psi
            if(isfield(TSE, 'VirtualCoilMartix'))
                if(~isempty(TSE.VirtualCoilMartix))
                    [TSE_sense_map, TSE.sense_Psi] = compress_sense_map_Psi(TSE.VirtualCoilMartix, TSE_sense_map,  TSE.sense_Psi);
                end
            end
            rs_command = sprintf('resize -c 0 %d', TSE.Ix_dim);
            TSE_sense_map = bart(rs_command, TSE_sense_map);
        end
        
        sense_map_3D = normalize_sense_map(TSE_sense_map(:,:,:,pars.enabled_ch ))+eps;
        sense_Psi = TSE.sense_Psi(pars.enabled_ch, pars.enabled_ch);
        
    else
        error('sense map source not indentified.')
    end
end

clear kspa_b0  %not useful afterwards

if(exist('sense_map_3D','var')&&exist('im_b0_ch_by_ch','var'))
    figure(21);
    slice = ceil(size(im_b0_ch_by_ch,3)/2);
    subplot(121); montage(abs(im_b0_ch_by_ch(:,:,slice,:)),'displayrange',[]); title('Check if they are match!'); xlabel('channel-by-channel');
    subplot(122); montage(abs(sense_map_3D(:,:,slice,:)),'displayrange',[]); xlabel('sense');
end

assert(sum(size(sense_map_3D)==[TSE.Ix_dim TSE.Iy_dim TSE.Iz_dim length(pars.enabled_ch)])==4,'SENSE map size is wrong!')

toc




%% preprocssing on phase error data
disp('Preprocessing on phase error data');
tic

%ref shot: when k0 being acquired
nonb0_idx = find(ismember(TSE.shot_matched, nonb0_shots_range));
ky_match_nonb0  = TSE.ky_matched(nonb0_idx);
kz_match_nonb0  = TSE.kz_matched(nonb0_idx);
shot_matched_nonb0 = TSE.shot_matched(nonb0_idx);
ref_shot_exact = unique(shot_matched_nonb0(find((ky_match_nonb0 == 0)&(kz_match_nonb0 == 0))));
ref_shot = find(ismember(ref_shot_exact, nonb0_shots_range))

clear  nonb0_idx ky_match_nonb0 kz_match_nonb0 shot_matched_nonb0 ref_shot_exact

%get nav_data;
nav_im_1 = double(nav_data(:,:,:,nonb0_shots_range)); %b0_shots_range for b0 correction; default: nonb0_shots_range

%smooth the "phase difference"
nav_im_1_diff = bsxfun(@rdivide, nav_im_1, nav_im_1(:,:,:,ref_shot)); %difference with the ref
nav_im_1_diff(find(isnan(nav_im_1_diff)))=0; nav_im_1_diff(find(isinf(nav_im_1_diff)))=0;

%set nav_im_1_diff phase to 0 for unreliable locations (optional and careful)
% unreliable_location = find(abs(nav_im_1)<1e4);
% nav_im_1_diff(unreliable_location)=abs(nav_im_1_diff(unreliable_location));


for sh = 1:size(nav_im_1_diff,4)
    nav_im_1_diff_phase_sm = smooth3(permute(angle(nav_im_1_diff(:,:,:,sh)),[3 1 2]),'box',pars.nav_phase_sm_kernel);
    nav_im_1_diff_phase_sm = permute(nav_im_1_diff_phase_sm, [2 3 1]);
    nav_im_1_diff_sm(:,:,:,sh) = abs(nav_im_1(:,:,:,sh)).*exp(1i.*nav_im_1_diff_phase_sm);
end
figure(501);
subplot(221); immontage4D(abs(nav_im_1),[]); colormap jet; title('mag map before sm');
subplot(223); immontage4D(angle(nav_im_1),[-pi pi]); colormap jet; title('phase map before sm');
subplot(222); immontage4D(angle(nav_im_1_diff),[-pi pi]); colormap jet; title('phase error before sm');
subplot(224); immontage4D(angle(nav_im_1_diff_sm),[-pi pi]); colormap jet; title('phase error after sm');
drawnow();

%interpolate to the correct size
kspace_interpo =  false;
if(kspace_interpo)
    assert(size(nav_im_1_diff_sm, 4)==nshot);
    
    nav_k_1 = bart('fft 7', nav_im_1_diff_sm);
    resize_command = sprintf('resize -c 0 %d 1 %d 2 %d 3 %d', TSE.Ix_dim/2, TSE.Ix_dim/2, TSE.Iz_dim, length(nonb0_shots_range)); %nav and TSE have the same FOV, but TSE have oversampling in x, so use ky instead of kx for the 1st dimension
    nav_k_1 = bart(resize_command, nav_k_1);
    nav_im_2 = bart('fft -i 7', nav_k_1);
    
    
else %linear intopolation
    [nav_x, ~, nav_z, nav_shots] = size(nav_im_1_diff_sm);
    nav_im_2 = zeros(TSE.Ix_dim, TSE.Ix_dim/2,  TSE.Iz_dim, length(nonb0_shots_range));
    padding_size = TSE.Ix_dim/2;
    padding_left = ceil(padding_size/2);
    x_loc_range = padding_left+1:padding_left+TSE.Ix_dim/2;
    
    if TSE.Iz_dim==1 %2D case; use imresize
        nav_im_2(x_loc_range,:,:,:) = imresize(nav_im_1_diff_sm, TSE.Iy_dim./size(nav_im_1_diff_sm,2));
    else % 3D use interp3
        
        for sh=1:nav_shots
            [X, Y, Z] = meshgrid(linspace(1,TSE.Ix_dim/2,nav_x),linspace(1,TSE.Ix_dim/2,nav_x),linspace(1,TSE.Iz_dim,nav_z) );  %corrdinate for original locations
            [Xq,Yq,Zq] = meshgrid(1:TSE.Ix_dim/2,1:TSE.Ix_dim/2,1:TSE.Iz_dim );  %corrdinate for intoplated locations
            nav_im_2(x_loc_range,:,:,sh) = interp3(X,Y,Z,squeeze(nav_im_1_diff_sm(:,:,:,sh)),Xq,Yq,Zq);
        end
        % nav_im_2 [TSE.Ix_dim, TSE.Ix_dim/2,  TSE.Iz_dim, n_shots]
        
        
    end
    
end


resize_command_2 = sprintf('resize -c 0 %d 1 %d 2 %d 3 %d', TSE.Ix_dim, TSE.Iy_dim, TSE.Iz_dim, length(nonb0_shots_range));
nav_im_2 = bart(resize_command_2, nav_im_2);
nav_im_2_masked = (bsxfun(@times, nav_im_2, TSE.sense_mask));  %mask


nav_im_3 = nav_im_2./abs(nav_im_2); %magnitude to 1;
nav_im_3(isnan(nav_im_3)) = 0; nav_im_3(isinf(nav_im_3)) = 0;
phase_error_3D = nav_im_3;

phase_error_3D = normalize_sense_map(phase_error_3D); %miss use normalize_sense_map
phase_error_3D = conj(phase_error_3D); %conj or not???
phase_error_3D = (bsxfun(@times, phase_error_3D, TSE.sense_mask));  %mask

% figure(5);
% immontage4D(angle(phase_error_3D),[-pi pi]); colormap jet; title('phase error maps int.')
% xlabel('shot #'); ylabel('slice locations');

assert(sum(size(phase_error_3D)==[TSE.Ix_dim TSE.Iy_dim TSE.Iz_dim length(nonb0_shots_range)])==4,'Phase Error map size is wrong!')

clear nav_im_3 nav_im_2_masked

toc


%% Preprocessing on kspace data for nonb0


disp('Preprocessing on kspace data for nonb0');
tic

% assert(sum(ismember(nonb0_shots_range, b0_shots_range))==0,'non_b0 shots overlap with b0 shots!')

% nonb0_shots_range = b0_shots_range; %backdoor: for calculating b0 images

% pad kx dim

padding_size = TSE.kx_dim-size(ima_k_spa_data,1);
padding_left = ceil(padding_size/2);
kx_sort_range = padding_left+1:padding_left+size(ima_k_spa_data,1);


ima_k_spa_data_kx_pad = zeros(TSE.kx_dim, size(ima_k_spa_data,2));
ima_k_spa_data_kx_pad(kx_sort_range,:)  = ima_k_spa_data;

ima_hybrid_k_spa_data = ifft1d(ima_k_spa_data_kx_pad);

clear ima_k_spa_data ima_k_spa_data_kx_pad

toc


%% msDWI recon
disp('recon');
image_corrected = zeros(TSE.Ix_dim, TSE.Iy_dim, TSE.Iz_dim);
tic;

% 3D recon (always 3D recon for large scale cases)


recon_x = pars.recon_x_locs;
%     recon_x = 300;
%
if(pars.parfor)        
    disp(['ParFor enabled; Othorgnalize kspa & sense not possible']);
    parfor x_idx = recon_x(1):recon_x(end)
        % for x_idx = recon_x(1):recon_x(end)
        
        disp(['current x location:',num2str(x_idx),'...']);
        %=========select phase_error ============================================
        phase_error = permute(permute(phase_error_3D(x_idx,:,:,:,:),[2 3 4 1]),[1 2 4 3]);
        %=========End ============================================
        
        if (sum(abs(phase_error(:))) == 0) % currect recon_x_loc is out of the SENSE make. nothing to recon
            image_corrected(x_idx,:,:) = zeros(TSE.Iy_dim, TSE.Iz_dim);
        else
            
            
            
            % [1]
            %=========select data. and sense map=====================================
            kspa = sort_k_spa_sh_by_sh_large_scale(ima_hybrid_k_spa_data(x_idx,:), nonb0_shots_range, TSE, pars);
            kspa = permute(kspa,[2 3 4 5 1]);
            sense_map = permute(sense_map_3D(x_idx,:,:,:),[2 3 4 1]);
            
            %========================================================================
            
            % [2]
            % ========Orthogonal SENSE maps: recombine coils to make the Psi map indentity mtx (SNR optimized)
            if(exist('sense_Psi', 'var'))
                [kspa, sense_map] = Ortho_data_sense(sense_Psi, kspa, sense_map);
                
            end
            % ========================================================================================
            
            
            % [3]
            image_corrected_current= msDWIrecon(kspa, sense_map, phase_error, pars.msDWIrecon);
            image_corrected_current(find(isnan(image_corrected_current))) = 0;
            image_corrected(x_idx,:,:) = image_corrected_current;
            
        end
        
    end
else %no parfor
    for x_idx = recon_x(1):recon_x(end)
        % for x_idx = recon_x(1):recon_x(end)
        
        disp(['current x location:',num2str(x_idx),'...']);
        %=========select phase_error ============================================
        phase_error = permute(permute(phase_error_3D(x_idx,:,:,:,:),[2 3 4 1]),[1 2 4 3]);
        %=========End ============================================
        
        if (sum(abs(phase_error(:))) == 0) % currect recon_x_loc is out of the SENSE make. nothing to recon
            image_corrected(x_idx,:,:) = zeros(TSE.Iy_dim, TSE.Iz_dim);
        else
            
            % display phase error------------------
            figure(411);
            montage(permute(squeeze(angle(phase_error)),[1 2 4 3]),'displayrange',[-pi pi]); colormap jet; title(['phase error map for x loc: ',num2str(x_idx)]);
            %---------------------------------------
            
            
            % [1]
            %=========select data. and sense map=====================================
            kspa = sort_k_spa_sh_by_sh_large_scale(ima_hybrid_k_spa_data(x_idx,:), nonb0_shots_range, TSE, pars);
            kspa = permute(kspa,[2 3 4 5 1]);
            sense_map = permute(sense_map_3D(x_idx,:,:,:),[2 3 4 1]);
            
            %========================================================================
            
            % [2]
            % ========Orthogonal SENSE maps: recombine coils to make the Psi map indentity mtx (SNR optimized)
                    if(exist('sense_Psi', 'var'))
                        [kspa, sense_map] = Ortho_data_sense(sense_Psi, kspa, sense_map);
            
                    end
            % ========================================================================================
            
            
            % [3]
            image_corrected(x_idx,:,:) = msDWIrecon(kspa, sense_map, phase_error, pars.msDWIrecon);
            image_corrected(x_idx,:,:) = 0;
            
        end
        
        %display
        figure(101);
        imshow(squeeze(abs(image_corrected(x_idx,:,:))),[]); title('msDWIrecon');  xlabel(['x loc: ',num2str(x_idx)]);
        drawnow();
        
        %check nav image & TSE image consistance
        figure(102);
        subplot(221);imshow(squeeze(abs(image_corrected(x_idx,:,:))),[]); title('reon');
        subplot(222);imshow(squeeze(mean(abs(nav_im_2(x_idx,:,:,:)),4)),[]); title('nav mean abs');
        subplot(223);imshow(squeeze(angle(image_corrected(x_idx,:,:))),[-pi pi]); title('angle recon');  xlabel(['x loc: ',num2str(x_idx)]);
        subplot(224);imshow(squeeze(angle(nav_im_2(x_idx,:,:,2))),[-pi pi]); title('nav angle');  xlabel(['x loc: ',num2str(x_idx)]);
        
        drawnow();
    end
end

%display

figure(109);
montage(permute(abs(image_corrected(80:250,:,:)),[1 2 4 3]),'displayrange',[]); title('msDWIrecon');


toc
end