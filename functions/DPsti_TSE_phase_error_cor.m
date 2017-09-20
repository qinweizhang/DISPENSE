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
%% Preprocessing on kspace data for b=0
disp('Preprocessing on kspace data for b=0');
if(isempty(pars.b0_shots))
    b0_shots_range = 1:shots_per_dyn; %by default, the first dynamic
else
    b0_shots_range = pars.b0_shots;
end

if(isempty(pars.nonb0_shots))
    nonb0_shots_range = 1: shots_per_dyn;
else
    nonb0_shots_range = pars.nonb0_shots;
end

assert(sum(ismember(nonb0_shots_range, b0_shots_range))==0,'non_b0 shots overlap with b0 shots!')


tic;
kx_dim = TSE.kxrange(2) - TSE.kxrange(1) + 1;
assert(kx_dim >=  size(ima_k_spa_data,1));
if(kx_dim>size(ima_k_spa_data,1))
    rs_command = sprintf('resize -c 0 %d ', kx_dim);
    ima_k_spa_data = bart(rs_command, ima_k_spa_data);
end

ky_dim = TSE.kyrange(2) - TSE.kyrange(1) + 1;  %max_ky * 2 + 1;
kz_dim = TSE.kzrange(2) - TSE.kzrange(1) + 1;
ch_dim = length(pars.enabled_ch); %TSE.ch_dim;
sh_dim = length([b0_shots_range nonb0_shots_range]); %range(TSE.shot_matched)+1;
kspa_all = zeros(kx_dim, ky_dim, kz_dim, ch_dim, sh_dim);

used_profile_nr = 0;
for prof_idx = 1:size(ima_k_spa_data, 2)
    ky_idx = TSE.ky_matched(prof_idx) + floor(ky_dim/2)+1;
    kz_idx = TSE.kz_matched(prof_idx) + floor(kz_dim/2)+1;
    ch_nr = mod(prof_idx, TSE.ch_dim) + 1;
    ch_idx = find(pars.enabled_ch==ch_nr);
    sh_nr = TSE.shot_matched(prof_idx);
    sh_idx = find([b0_shots_range nonb0_shots_range]==sh_nr);
    if(~isempty(ch_idx)&&~isempty(sh_idx))
        used_profile_nr = used_profile_nr+1;
        kspa_all(:,ky_idx,kz_idx,ch_idx, sh_idx) = ima_k_spa_data(:,prof_idx);
    end
end

temp = kspa_all(round(kx_dim/2),:,:,:,:);
assert(used_profile_nr == sum(abs(temp(:))>0)); clear temp
size(kspa_all)


% remove stupid checkerboard pattern
che=create_checkerboard([1,size(kspa_all,2),size(kspa_all,3)]);
kspa_all=bsxfun(@times,kspa_all,che);


kspa_b0 = sum(kspa_all(:,:,:,:,[1:length(b0_shots_range)]), 5); %4D b0 kspace [kx ky kz nc]

pp=bart('fft -i 7',kspa_b0);
im_b0=bart('rss 8',pp);
figure(1); montage(permute(abs(im_b0),[1 2 4 3]),'displayrange',[]); title('b0 images');
clear pp
toc;

%% Preprocessing on sense data
disp('Preprocessing on sense data');
tic

% %sense mask || now it should be calculated outside and stored in TSE structure
if(isfield(TSE, 'sense_mask'))
    if(isempty(TSE.sense_mask)) %if empty calc again
        dim = [size(kspa_b0, 2) size(kspa_b0, 2) size(kspa_b0,3)];
        os = [1, 1, 1];
        sense_map_temp = get_sense_map_external(pars.sense_ref, pars.data_fn, pars.coil_survey, dim, os);
        rs_command = sprintf('resize -c 0 %d', size(kspa_b0, 1));
        sense_map_temp = bart(rs_command, sense_map_temp);
        
        TSE.sense_mask = abs(sense_map_temp(:,:,:,1 ))>0;
        clear sense_map_temp;
    end
else %if not exist, calc again
    dim = [size(kspa_b0, 2) size(kspa_b0, 2) size(kspa_b0,3)];
    os = [1, 1, 1];
    sense_map_temp = get_sense_map_external(pars.sense_ref, pars.data_fn, pars.coil_survey, dim, os);
    rs_command = sprintf('resize -c 0 %d', size(kspa_b0, 1));
    sense_map_temp = bart(rs_command, sense_map_temp);
    
    TSE.sense_mask = abs(sense_map_temp(:,:,:,1 ))>0;
    clear sense_map_temp;
end

%sense maps
if(length(pars.enabled_ch)==1) %one channel
    warning('Kerry: This is one channel recon!')
    sense_map_3D = ones(size(kspa_all,1),size(kspa_all,2),size(kspa_all,3));
        
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
        
        
        %match the size of TSE_sens_map to kspa_xyz
        %TODO
    elseif(strcmp(pars.sense_map, 'external'))
        if(isempty(TSE_sense_map))
            dim = [size(kspa_b0, 2) size(kspa_b0, 2) size(kspa_b0,3)];
            os = [1, 1, 1];
            TSE_sense_map = get_sense_map_external(pars.sense_ref, pars.data_fn, pars.coil_survey, dim, os);
            rs_command = sprintf('resize -c 0 %d', size(kspa_b0, 1));
            TSE_sense_map = bart(rs_command, TSE_sense_map);
        end
        sense_map_3D = normalize_sense_map(TSE_sense_map(:,:,:,pars.enabled_ch ))+eps;
    else
        error('sense map source not indentified.')
    end
end
toc

%% Preprocessing on kspace data for nonb0
disp('Preprocessing on kspace data for nonb0');

tic
kspa_xyz = kspa_all(:,:,:,:,(length(b0_shots_range)+1):end);

kk = sum(kspa_xyz, 5); %4D b0 kspace [kx ky kz nc]
pp=bart('fft -i 7',kk);
im_nonb0=bart('rss 8',pp);
figure(2); montage(permute(abs(im_nonb0),[1 2 4 3]),'displayrange',[]); title('direct recon');
clear kk pp

%to hybrid space
clear kspa_x_yz;
kspa_x_yz = ifft1d(kspa_xyz);
toc



%% preprocssing on phase error data
disp('Preprocessing on phase error data');
tic
%match size of nav_im and to kspa_xyz
nav_im_1 = double(nav_data(:,:,:,nonb0_shots_range));
[kx, ky, kz, nc, nshot] = size(kspa_xyz);
assert(size(nav_im_1, 4)==nshot);

nav_k_1 = bart('fft 7', nav_im_1);
resize_command = sprintf('resize -c 0 %d 1 %d 2 %d 3 %d', ky, ky, kz, nshot); %nav and TSE have the same FOV, but TSE have oversampling in x, so use ky instead of kx for the 1st dimension
nav_k_1 = bart(resize_command, nav_k_1);
nav_im_2 = bart('fft -i 7', nav_k_1);

resize_command_2 = sprintf('resize -c 0 %d 1 %d 2 %d 3 %d', kx, ky, kz, nshot);
nav_im_2 = bart(resize_command_2, nav_im_2);
nav_im_2 = nav_im_2./abs(nav_im_2); %magnitude to 1;
nav_im_2(isnan(nav_im_2)) = 0; nav_im_2(isinf(nav_im_2)) = 0;

phase_error_3D = nav_im_2;
k0_idx = [floor(kx_dim/2)+1 floor(ky_dim/2)+1 floor(kz_dim/2)+1];

for sh =1:nshot
    if(abs(kspa_xyz(k0_idx(1),k0_idx(2),k0_idx(3), 1, sh))>0)
        ref_shot = sh;
    end
end
phase_error_3D = bsxfun(@rdivide, phase_error_3D, phase_error_3D(:,:,:,ref_shot)); %difference with the ref
phase_error_3D = normalize_sense_map(phase_error_3D); %miss use normalize_sense_map
% phase_error_3D = conj(bsxfun(@times, phase_error_3D, abs(sense_map_3D(:,:,:,1))>eps)) + eps; %mask the outer region
phase_error_3D = conj(bsxfun(@times, phase_error_3D, TSE.sense_mask));  %conj or not???

figure(5);
immontage4D(angle(phase_error_3D),[-pi pi]); colormap jet; title('phase error maps int.')
xlabel('shot #'); ylabel('slice locations');


toc
%% msDWI recon
disp('recon');
image_corrected = zeros(kx, ky, kz);
tic;
if (kz>1)
    %% 3D recon
    recon_x = [100: 220]; %pars.recon_x_locs;
    for x_idx = 1:length(recon_x)
        
        recon_x_loc = recon_x(x_idx);
        
        %=========select data. fixed=============================================
        kspa = permute(kspa_x_yz(recon_x_loc, :, :, :, :),[2 3 4 5 1]);
        sense_map = permute(sense_map_3D(recon_x_loc,:,:,:),[2 3 4 1]);
        phase_error = permute(permute(phase_error_3D(recon_x_loc,:,:,:,:),[2 3 4 1]),[1 2 4 3]);
        %========================================================================
        
        image_corrected(recon_x_loc,:,:) = msDWIrecon(kspa, sense_map, phase_error, pars.msDWIrecon);
        
        image_corrected(isnan(image_corrected)) = 0;
        
        %display
        figure(101);
        subplot(131);imshow(squeeze(abs(im_b0(recon_x_loc,:,:))),[]); title('b0');
        subplot(132);imshow(squeeze(abs(im_nonb0(recon_x_loc,:,:))),[]); title('direct recon');
        subplot(133);imshow(squeeze(abs(image_corrected(recon_x_loc,:,:))),[]); title('msDWIrecon');
        
    end
    
    toc;
    
    %display
    figure(109);
    subplot(141); montage(permute(abs(im_b0(:,:,:)),[1 2 4 3]),'displayrange',[]); title('b0');
    subplot(142); montage(permute(abs(im_nonb0(:,:,:)),[1 2 4 3]),'displayrange',[]); title('direct recon');
    subplot(143); montage(permute(abs(image_corrected(:,:,:)),[1 2 4 3]),'displayrange',[]); title('msDWIrecon');
else
    %% 2D recon: remove 3rd dimension
    
        image_corrected = msDWIrecon(permute(kspa_xyz,[1 2 4 5 3]), squeeze(sense_map_3D), phase_error_3D, pars.msDWIrecon);
%     image_corrected = msDWIrecon(permute(kspa_all(:,:,:,:,[1:length(b0_shots_range)]),[1 2 4 5 3]), squeeze(sense_map_3D), phase_error_3D, pars.msDWIrecon);
    %display
    figure(110);
    subplot(141); imshow(abs(im_b0),[]); title('b0');
    subplot(142);  imshow(abs(im_nonb0),[]); title('b0'); title('direct recon');
    subplot(144);  imshow(abs(image_corrected),[]); title('b0'); title('msDWIrecon');
    
    
    
end

end