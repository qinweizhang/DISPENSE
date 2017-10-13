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
function image_corrected = DPsti_TSE_phase_error_cor(ima_k_spa_data, TSE, TSE_sense_map, nav_data, pars)

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

%% Preprocessing on kspace data for b=0

disp('Preprocessing on kspace data for b=0');
if(isempty(pars.b0_shots))
    b0_shots_range = 1:shots_per_dyn; %by default, the first dynamic
else
    b0_shots_range = pars.b0_shots;
end

tic;

kspa_b0 = sort_k_spa_sh_by_sh(ima_k_spa_data, b0_shots_range, TSE, pars);

kspa_b0_combshot = sum(kspa_b0, 5)./sum(abs(kspa_b0)>0, 5); %4D b0 kspace [kx ky kz nc]; non-zeros average
kspa_b0_combshot(find(isnan(kspa_b0_combshot)))=0; kspa_b0_combshot(find(isinf(kspa_b0_combshot)))=0;

im_b0_ch_by_ch=ifft3d(kspa_b0_combshot);
im_b0=sqrt(sum(abs(im_b0_ch_by_ch).^2, 4));
figure(1); montage(permute(abs(im_b0),[1 2 4 3]),'displayrange',[]); title('b0 images');

clear kspa_b0_combshot

toc;

%% Preprocessing on kspace data for nonb0

disp('Preprocessing on kspace data for nonb0');

if(isempty(pars.nonb0_shots))
    nonb0_shots_range = 1: shots_per_dyn;
else
    nonb0_shots_range = pars.nonb0_shots;
end

% assert(sum(ismember(nonb0_shots_range, b0_shots_range))==0,'non_b0 shots overlap with b0 shots!')

% nonb0_shots_range = b0_shots_range; %backdoor: for calculating b0 images

tic;

kspa_xyz = sort_k_spa_sh_by_sh(ima_k_spa_data, nonb0_shots_range, TSE, pars);

kk = sum(kspa_xyz, 5)./ sum(abs(kspa_xyz)>0, 5); %4D b0 kspace [kx ky kz nc]; non-zero average
kk(find(isnan(kk)))=0; kk(find(isinf(kk)))=0;

im_nonb0_ch_by_ch=ifft3d(kk);
im_nonb0=sqrt(sum(abs(im_nonb0_ch_by_ch).^2, 4));
figure(2); montage(permute(abs(im_nonb0),[1 2 4 3]),'displayrange',[]); title('direct recon');
clear kk pp

%to hybrid space for 3D cases
if(TSE.kz_dim>1)
    clear kspa_x_yz;
    kspa_x_yz = ifft1d(kspa_xyz);
    clear kspa_xyz
end

toc

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
k0_idx = [floor(TSE.kx_dim/2)+1 floor(TSE.ky_dim/2)+1 floor(TSE.kz_dim/2)+1];
for sh =1:length(nonb0_shots_range)
    if(exist('kspa_xyz','var'))
        if(abs(kspa_xyz(k0_idx(1),k0_idx(2),k0_idx(3), 1, sh))>0)
            ref_shot = sh;
        end
    else
        if(abs(kspa_x_yz(k0_idx(1),k0_idx(2),k0_idx(3), 1, sh))>0)
            ref_shot = sh;
        end
    end
    
end

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


nav_im_2 = nav_im_2./abs(nav_im_2); %magnitude to 1;
nav_im_2(isnan(nav_im_2)) = 0; nav_im_2(isinf(nav_im_2)) = 0;
phase_error_3D = nav_im_2;

phase_error_3D = normalize_sense_map(phase_error_3D); %miss use normalize_sense_map
phase_error_3D = conj(phase_error_3D); %conj or not???
phase_error_3D = (bsxfun(@times, phase_error_3D, TSE.sense_mask));  %mask

figure(5);
immontage4D(angle(phase_error_3D),[-pi pi]); colormap jet; title('phase error maps int.')
xlabel('shot #'); ylabel('slice locations');

assert(sum(size(phase_error_3D)==[TSE.Ix_dim TSE.Iy_dim TSE.Iz_dim length(nonb0_shots_range)])==4,'Phase Error map size is wrong!')

toc
%% msDWI recon
disp('recon');
image_corrected = zeros(TSE.Ix_dim, TSE.Iy_dim, TSE.Iz_dim);
tic;
if (TSE.Iz_dim>1)
    %% 3D recon
    recon_x = pars.recon_x_locs;
%     recon_x = 150;
    for x_idx = 1:length(recon_x)
        
        recon_x_loc = recon_x(x_idx);
        
        %=========select data. fixed=============================================
        kspa = permute(kspa_x_yz(recon_x_loc, :, :, :, :),[2 3 4 5 1]);
        sense_map = permute(sense_map_3D(recon_x_loc,:,:,:),[2 3 4 1]);
        phase_error = permute(permute(phase_error_3D(recon_x_loc,:,:,:,:),[2 3 4 1]),[1 2 4 3]);
        %========================================================================
        
        % ========Orthogonal SENSE maps: recombine coils to make the Psi map indentity mtx (SNR optimized)
        if(exist('sense_Psi', 'var'))
            L = chol(sense_Psi,'lower'); %Cholesky decomposition; L: lower triangle
            L_inv = inv(L);
            for c = 1:size(sense_Psi,1)
                %recombine sense map
                sense_map_orthocoil(:,:,c) = sum(bsxfun(@times, sense_map, permute(L_inv(c,:),[1 3 2])), 3);
                %recombine kspa map
                kspa_orthocoil(:,:,c,:) = sum(bsxfun(@times, kspa, permute(L_inv(c,:),[1 3 2])), 3);
            end
            figure(401);
            subplot(121); montage(permute(abs(sense_map),[1 2 4 3]),'displayrange',[]); title('originial SENSE map')
            subplot(122); montage(permute(abs(sense_map_orthocoil),[1 2 4 3]),'displayrange',[]); title('Orthogonal SENSE map')
            
            sense_map =  sense_map_orthocoil;
            kspa = kspa_orthocoil;
            clear sense_map_orthocoil kspa_orthocoil
            %renormalize sense
            sense_map = squeeze(normalize_sense_map(permute(sense_map,[1 2 4 3])));
        end
        % ========================================================================================
        
        
        image_corrected(recon_x_loc,:,:) = msDWIrecon(kspa, sense_map, phase_error, pars.msDWIrecon);
        
        image_corrected(isnan(image_corrected)) = 0;
        
        %display
        figure(101);
        if(exist('im_b0'))
            subplot(131);imshow(squeeze(abs(im_b0(recon_x_loc,:,:))),[]); title('b0');
        end
        subplot(132);imshow(squeeze(abs(im_nonb0(recon_x_loc,:,:))),[]); title('direct recon');
        subplot(133);imshow(squeeze(abs(image_corrected(recon_x_loc,:,:))),[]); title('msDWIrecon');  xlabel(['x loc: ',num2str(recon_x_loc)]);
        drawnow();
        %         figure(102); montage(angle((phase_error)),[-pi pi]); colormap jet
        
        
    end
    
    
    %display
    if(1) %big screen
        figure(109);
        if(exist('im_b0','var'))
            subplot(141); montage(permute(abs(im_b0(80:250,:,:)),[1 2 4 3]),'displayrange',[]); title('b0');
        end
        subplot(142); montage(permute(abs(im_nonb0(80:250,:,:)),[1 2 4 3]),'displayrange',[]); title('direct recon');
        subplot(143); montage(permute(abs(image_corrected(80:250,:,:)),[1 2 4 3]),'displayrange',[]); title('msDWIrecon');
    else
        if(exist('im_b0','var'))
            figure(109); montage(permute(abs(im_b0(80:250,:,:)),[1 2 4 3]),'displayrange',[]); title('b0');
        end
        figure(110); montage(permute(abs(im_nonb0(80:250,:,:)),[1 2 4 3]),'displayrange',[]); title('direct recon');
        figure(111); montage(permute(abs(image_corrected(80:250,:,:)),[1 2 4 3]),'displayrange',[]); title('msDWIrecon');
    end
else
    %% 2D recon: remove 3rd dimension
    
    
    ksap = permute(kspa_xyz,[1 2 4 5 3]);
    sense_map = squeeze(sense_map_3D);
    
    % ========Orthogonal SENSE maps: recombine coils to make the Psi map indentity mtx (SNR optimized)
    if(exist('sense_Psi', 'var'))
        L = chol(sense_Psi,'lower'); %Cholesky decomposition; L: lower triangle
        L_inv = inv(L);
        for c = 1:size(sense_Psi,1)
            %recombine sense map
            sense_map_orthocoil(:,:,c) = sum(bsxfun(@times, sense_map, permute(L_inv(c,:),[1 3 2])), 3);
            %recombine kspa map
            kspa_orthocoil(:,:,c,:) = sum(bsxfun(@times, kspa, permute(L_inv(c,:),[1 3 2])), 3);
        end
        figure(401);
        subplot(121); montage(permute(abs(sense_map),[1 2 4 3]),'displayrange',[]); title('originial SENSE map')
        subplot(122); montage(permute(abs(sense_map_orthocoil),[1 2 4 3]),'displayrange',[]); title('Orthogonal SENSE map')
        
        sense_map =  sense_map_orthocoil;
        kspa = kspa_orthocoil;
        clear sense_map_orthocoil kspa_orthocoil
        %renormalize sense
        sense_map = squeeze(normalize_sense_map(permute(sense_map,[1 2 4 3])));
    end
    % ========================================================================================
    
    
    image_corrected = msDWIrecon(ksap, sense_map, phase_error_3D, pars.msDWIrecon);
    %     image_corrected = msDWIrecon(permute(kspa_all(:,:,:,:,[1:length(b0_shots_range)]),[1 2 4 5 3]), squeeze(sense_map_3D), phase_error_3D, pars.msDWIrecon);
    %display
    figure(120);
    subplot(141); imshow(abs(im_b0),[]); title('b0');
    subplot(142);  imshow(abs(im_nonb0),[]); title('b0'); title('direct recon');
    subplot(143);  imshow(abs(image_corrected),[]); title('b0'); title('msDWIrecon');
    
    
    
end
toc
end