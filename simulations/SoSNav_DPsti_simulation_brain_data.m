%==========DPsti-TSE nonrigid motion correction simulation============
cd('/home/qzhang/lood_storage/divi/Users/qzhang/SoSND/SoSND/simulations/brain_test_data')
clear; close all; clc
load('test_data_09_24.mat')  %  test_data_09_24



tic;
noiselevel = 0;


%% STEP1: generate sampling pattern

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
LRT_nav_mask_enable = true;


load_default = false;
load_sparse_profile = ~load_default;


if(load_default)
    %% low-high radial sampling pattern; from the TSE data
    sh_dim = max(TSE.shot_matched) - min(TSE.shot_matched) + 1;
    
    ky_labels = TSE.ky_matched + ceil(TSE.ky_dim/2); %0-->TSE.ky_dim/2
    ky_labels = ky_labels(1:TSE.ch_dim:end);
    kz_labels = TSE.kz_matched +  ceil(TSE.kz_dim/2); %0-->TSE.kz_dim/2
    kz_labels = kz_labels(1:TSE.ch_dim:end);
    prof_per_shot = length(ky_labels)/sh_dim;
    
elseif(load_sparse_profile)
    %% sparse sampling pattern
    sh_dim = 42;
    ky_labels = sparse_ky_profile + ceil(TSE.ky_dim/2); %0-->TSE.ky_dim/2
    kz_labels = sparse_kz_profile + ceil(TSE.kz_dim/2); %0-->TSE.ky_dim/2
    prof_per_shot = length(ky_labels)/sh_dim;
end

mask = [];

for sh = 1:sh_dim
    range = (sh-1)*prof_per_shot+1:sh*prof_per_shot;
    m = sparse(ky_labels(range), kz_labels(range), ones(prof_per_shot,1), TSE.ky_dim, TSE.kz_dim);
    mask(1,:,:,1,sh) =full(m);
end
ttt(:,:,1,:) = squeeze(mask);
disp([num2str(sum(mask(:))), ' profiles in total.']);

k0_idx = [floor(TSE.kx_dim/2)+1 ceil(TSE.ky_dim/2) ceil(TSE.kz_dim/2)];

%plot
figure(12);
subplot(141); montage(ttt, 'displayrange',[]);xlabel('kz'); ylabel('ky'); title('sampling mask shot-by-shot');
subplot(142); imagesc(sum(ttt,4), [0 2]); colormap jet; axis off; axis equal; xlabel('kz'); ylabel('ky'); title('sample times'); colorbar;
ttt = sum(bsxfun(@times, ttt, permute(1:42, [4 3 1 2])),4);
subplot(143); imagesc(ttt, [0 sh_dim]); colormap jet; axis off; axis equal; xlabel('kz'); ylabel('ky'); title('smaple shot #'); colorbar;
clear ttt;

%================separate navigator for LRT============================
LRT_nav_mask = [];
if(LRT_nav_mask_enable)
    LRT_nav_mask = zeros(size(mask));
    LRT_nav_mask(:,k0_idx(2)-5:k0_idx(2)+5,k0_idx(3)-5:k0_idx(3)+5,:,:) = 1;
end

%%  STEP2: generate phase error for every shot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




generate_phase_error = false;
load_phase_error = ~generate_phase_error;


if(generate_phase_error)
    phase_error_3D = [];
    for shot = 1:sh_dim
        
        %% global phase error
        %         pe_temp = exp(1i .* shot .* ones(TSE.kx_dim, TSE.ky_dim, TSE.kz_dim));
        %         phase_error_3D = cat(5, phase_error_3D, pe_temp);
        %% linear phase error
        %         pe_2D = linear2Dmap(shot*3, 0.2, TSE.kx_dim, TSE.ky_dim);
        %         pe_temp = exp(1i .* pe_2D);
        %         phase_error_3D = cat(5, phase_error_3D, pe_temp);
        
        %% random phase error
        rand_phase = 50 * pi *random_phase_map(TSE.kx_dim, TSE.ky_dim, TSE.kz_dim,1);
        pe_temp = exp(1i.*rand_phase);
        phase_error_3D = cat(5, phase_error_3D, pe_temp);
        
    end
    if(size(phase_error_3D, 3)~=TSE.kz_dim)
        phase_error_3D = repmat(phase_error_3D, [1 1 TSE.kz_dim 1 1 ]);
    end
elseif(load_phase_error)
    %% external phaes error
    assert(sum(size(phase_error_3D_stored) == [TSE.kx_dim, TSE.ky_dim, TSE.kz_dim, sh_dim]) == 4 );
    phase_error_3D = permute(phase_error_3D_stored, [1 2 3 5 4]);
    size(phase_error_3D)
end

phase_error_3D  = phase_error_3D./abs(phase_error_3D);
phase_error_3D(find(isnan(phase_error_3D))) = 0;
phase_error_3D(find(isinf(phase_error_3D))) = 0;
















%% STEP3: corrupted im_ref + STEP4: recon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_corrected = zeros(size(im_ref));
im_recon_direct = zeros(size(im_ref));

ch_id = [3 5  6 10 15 20];
recon_shot_range = 1:20;
recon_x_range = 150; %130:170;
for x = 1:length(recon_x_range)
    
    recon_x_loc = recon_x_range(x);
    
    %recon parameter------------------------------------------------
    pars = initial_msDWIrecon_Pars;
    
    pars.method='LRT'; %POCS_ICE CG_SENSE_I CG_SENSE_K LRT
    
    pars.CG_SENSE_I.lamda=2;
    pars.CG_SENSE_I.nit=5;
    pars.CG_SENSE_I.tol = 1e-10;
    
    pars.POCS.Wsize = [10 10];  %no point to be bigger than navigator area
    pars.POCS.nit = 300;
    pars.POCS.tol = 1e-10;
    pars.POCS.lamda = 1;
    pars.POCS.nufft = false;
    
    pars.LRT.Lg=1;
    pars.LRT.L3=6;
    pars.LRT.L4=2;
    pars.LRT.mu = 1.5;
    pars.LRT.beta = 1;
    pars.LRT.lambda = 3e-2;
    
    pars.LRT.sparsity_transform='TV';
    pars.LRT.Imref=repmat(squeeze(im_ref(recon_x_loc,:,:,:)), [1 1 length(recon_shot_range)+1 2]);
    pars.LRT.x=15;
    pars.LRT.y=88;
    pars.LRT.increase_penalty_parameters=false;
    pars.LRT.inspectLg=false;
    pars.LRT.subspacedim1=1;
    pars.LRT.subspacedim2=1;
    pars.LRT.G.precon=true;
    pars.LRT.G.maxiter = 10;
    pars.LRT.scaleksp=1;
    pars.LRT.niter = 20;
    
    
    %--------------------------------------------------------------
    
    %% corrupt data
    
    %select channel ids

    if(length(ch_id)>1)
        sense_map_3D = normalize_sense_map(TSE_sense_map(recon_x_loc,:,:,ch_id));
    else
        sense_map_3D = ones(size(im_ref));
    end
    
    %select pe x loc
    phase_error_3D_current = phase_error_3D(recon_x_loc,:,:,:,:);
    
    
    
    im_ref_ch =  bsxfun(@times, im_ref(recon_x_loc,:,:), sense_map_3D);
    s = ceil(TSE.kz_dim/2);
    %     figure(121); montage(permute(squeeze(abs(im_ref_ch)),[1 2 4 3]),'displayrange',[]); title('reference image, channel by channel');
    im_pe = bsxfun(@times,im_ref_ch,phase_error_3D_current);
    %======================================================================
    
    %================self_navigator size for POCS==========================
    if(strcmp( pars.method, 'POCS_ICE')) % for LRT separate navigator
        
        mask(:,k0_idx(2)-5:k0_idx(2)+5,k0_idx(3)-5:k0_idx(3)+5,:,:) = 1;
    end
    
    
    %================SENSE undersampling===================================
%         mask(:,1:2:end,:,:,:) = 0;
    
    kspa_xyz = bsxfun(@times, fft3d(im_pe), mask);
    if(~isempty(LRT_nav_mask)) % for LRT separate navigator
        kspa_xyz_nav =  bsxfun(@times, fft3d(im_pe), LRT_nav_mask);
    end
    %     kspa_xyz = fft3d(im_pe);
    clear im_pe;
    
    %add noise
    
    if noiselevel>0
        kspa_xyz=kspa_xyz+(randn(size(kspa_xyz)).*mean(kspa_xyz(find(abs(kspa_xyz)>0))).*noiselevel).*(abs(kspa_xyz)>0);
        if(~isempty(LRT_nav_mask)) % for LRT separate navigator
            kspa_xyz_nav = kspa_xyz_nav+(randn(size(kspa_xyz_nav)).*mean(kspa_xyz_nav(find(abs(kspa_xyz_nav)>0))).*noiselevel).*(abs(kspa_xyz_nav)>0);
        end
    end
    
    
    %to hybrid space
    clear kspa_x_yz;
    kspa_x_yz = ifft1d(kspa_xyz);
    if(~isempty(LRT_nav_mask)) % for LRT separate navigator
        kspa_x_yz_nav = ifft1d(kspa_xyz_nav);
    end
    
    
    %direct recon
    k_combine_shots=sum(kspa_xyz,5)./sum(abs(kspa_xyz)>0,5); %non zeros average
    k_combine_shots(find(isnan(k_combine_shots)))=0; k_combine_shots(find(isinf(k_combine_shots)))=0;
    
    im_ch_by_ch=ifft3d(k_combine_shots);
    im_recon_direct(recon_x_loc, :,:)=sqrt(sum(abs(im_ch_by_ch).^2, 4));
    %     figure(4); montage(permute(abs(im_recon_direct),[1 2 4 3]),'displayrange',[]); title('direct recon');
    clear  im_ch_by_ch k_combine_shots
    
    if(~isempty(LRT_nav_mask)) % for LRT separate navigator
        kk=sum(kspa_xyz_nav,5)./sum(abs(kspa_xyz_nav)>0,5); %non zeros average
        kk(find(isnan(kk)))=0; kk(find(isinf(kk)))=0;
        
        kk_c=ifft3d(kk);
        pp=sqrt(sum(abs(kk_c).^2, 4));
        figure(41); montage(permute(abs(pp),[2 3 4 1]),'displayrange',[]); title('direct of navigator');
        clear  pp ll kk_c
    end
    
    
    toc;
    
    %% recon
    
    tic;
    
    
    % preprocss phase_error_maps for CG_CENSE_I
    for sh =1:sh_dim
        if(abs(kspa_xyz(ceil(size(kspa_xyz,1)/2),k0_idx(2),k0_idx(3), 1, sh))>0)
            ref_shot = sh;
        end
    end
    
    
    phase_error_3D_current = squeeze(bsxfun(@rdivide, phase_error_3D_current, phase_error_3D_current(:,:,:,ref_shot))); %difference with the first
    phase_error_3D_current(find(isnan(phase_error_3D_current)))=0;     phase_error_3D_current(find(isinf(phase_error_3D_current)))=0;
    phase_error_3D_current = permute(normalize_sense_map(permute(phase_error_3D_current,[1 2 4 3])),[5 1 2 3 4]); %miss use normalize_sense_map
    
    
    %=========select data. fixed=============================================
    kspa = permute(kspa_x_yz(1, :, :, :, recon_shot_range),[2 3 4 5 1]);
    kspa_bart =  sum(kspa,4)./sum(abs(kspa)>0,4); kspa_bart(find(isnan(kspa_bart)))=0; kspa_bart(find(isinf(kspa_bart)))=0;
    if(strcmp( pars.method, 'LRT')) % for LRT separate navigator
        %cat navigator 
        kspa = cat(5, permute(kspa_x_yz_nav(1, :, :, :,  recon_shot_range),[2 3 4 5 1]), kspa);
        %cat "b0" info into the first row
        kspa_ref_ch = fft2d(squeeze(im_ref_ch))*10; %b0 image has higher intensity
        nav_ref_ch = kspa_ref_ch.*squeeze(abs(kspa_x_yz_nav(:,:,:,:,1))>0);
        kspa = cat(4, cat(5, nav_ref_ch, kspa_ref_ch), kspa);
        
    end
    sense_map = permute(sense_map_3D(1,:,:,:),[2 3 4 1])+eps;
    phase_error = permute(squeeze(phase_error_3D_current(1,:,:,:, recon_shot_range)),[1 2 4 3])+eps;
    %========================================================================
    
    
    
    image_corrected(recon_x_loc,:,:) = msDWIrecon(kspa, (sense_map), (phase_error), pars);
    
    
    %display
    
    figure(102);
    subplot(141);imshow(squeeze(abs(im_ref(recon_x_loc,:,:))),[]); title('reference');
    subplot(142);imshow(squeeze(abs(im_recon_direct(recon_x_loc,:,:))),[]); title('direct recon');
    subplot(143);imshow(squeeze(abs(image_corrected(recon_x_loc,:,:))),[]); title('msDWIrecon'); xlabel(pars.method);
    if(strcmp( pars.method, 'LRT')) 
        bart_recon(recon_x_loc,:,:) = bart('pics -l2 -r0.00001 -i200', permute(kspa_bart, [1 2 4 3]), permute(sense_map,[1 2 4 3]));
        subplot(144);imshow(squeeze(abs(bart_recon(recon_x_loc,:,:))),[]); title('bart PICS');
    end
    drawnow();
end
figure(109);
% subplot(141);
montage(permute(abs(im_ref(recon_x_range,:,:)),[1 2 4 3]),'displayrange',[]); title('reference');
figure(110);
%     subplot(142);
montage(permute(abs(im_recon_direct(recon_x_range,:,:)),[1 2 4 3]),'displayrange',[]); title('direct recon');
figure(111);
%     subplot(143);
montage(permute(abs(image_corrected(recon_x_range,:,:)),[1 2 4 3]),'displayrange',[]); title('msDWIrecon');

toc;