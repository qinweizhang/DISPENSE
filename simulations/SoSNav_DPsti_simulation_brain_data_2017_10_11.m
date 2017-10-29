%==========DPsti-TSE nonrigid motion correction simulation============
cd('/home/qzhang/lood_storage/divi/Users/qzhang/SoSND/SoSND/simulations/brain_test_data')
clear; close all; clc
load('test_data_10_11.mat')  %  test_data_09_24

%%  STEP1: generate phase error for every shot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




generate_phase_error = false;
load_phase_error = ~generate_phase_error;



sh_dim = max(TSE.shot_matched) - min(TSE.shot_matched) + 1;
if(generate_phase_error)
    phase_error_3D = [];
    for shot = 1:sh_dim
        
        %% global phase error
        %                 pe_temp = exp(1i .* shot .* ones(TSE.kx_dim, TSE.ky_dim, TSE.kz_dim));
        %                 phase_error_3D = cat(5, phase_error_3D, pe_temp);
        %% linear phase error
        %                 pe_2D = linear2Dmap(shot*3, 0.2, TSE.kx_dim, TSE.ky_dim);
        %                 pe_temp = exp(1i .* pe_2D);
        %                 phase_error_3D = cat(5, phase_error_3D, pe_temp);
        
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

phase_error_3D  = phase_error_3D./abs(phase_error_3D); %magnitude=1
phase_error_3D = phase_error_3D.*exp(1i*angle(phase_error_3D*5)); %more phase error
phase_error_3D(find(isnan(phase_error_3D))) = 0;
phase_error_3D(find(isinf(phase_error_3D))) = 0;
















%% STEP2: update mask + STEP3: corrupted im_ref + STEP4: recon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all



load_default = false;
load_sparse_profile = ~load_default;

noiselevel = 0;


image_corrected = zeros(size(im_ref));
im_recon_direct = zeros(size(im_ref));
simulate_ky_sense_2 = false;
simulate_kz_sense_2 = false;

ch_id = [1:10];
recon_shot_range = 1:26;
recon_x_range = 190; %130:170;
for x = 1:length(recon_x_range)
    %% Set recon parameters
    
    recon_x_loc = recon_x_range(x);
    
    %recon parameter------------------------------------------------
    pars = initial_msDWIrecon_Pars;
    
    pars.method='LRT'; %POCS_ICE CG_SENSE_I CG_SENSE_K LRT
    
    pars.CG_SENSE_I.lamda=0.1;
    pars.CG_SENSE_I.nit=30;
    pars.CG_SENSE_I.tol = 1e-40;
    
    pars.POCS.Wsize = [10 10];  %no point to be bigger than navigator area
    pars.POCS.nit = 300;
    pars.POCS.tol = 1e-10;
    pars.POCS.lamda = 0;
    pars.POCS.nufft = false;
    
    pars.LRT.Lg=3;
    pars.LRT.L3=4;
    pars.LRT.L4=1;
    pars.LRT.mu = 2e2;
    pars.LRT.beta = 1;
    pars.LRT.lambda = 2e-2;
    
    pars.LRT.sparsity_transform='TV';
    pars.LRT.Imref=cat(3, repmat(squeeze(im_b0_ref(recon_x_loc,:,:,:)), [1 1 1 2]), repmat(squeeze(im_ref(recon_x_loc,:,:,:)), [1 1 length(recon_shot_range) 2]));
    pars.LRT.x=15;
    pars.LRT.y=88;
    pars.LRT.increase_penalty_parameters=false;
    pars.LRT.inspectLg=false;
    pars.LRT.subspacedim1=1;
    pars.LRT.subspacedim2=1;
    pars.LRT.G.precon=true;
    pars.LRT.G.maxiter = 10;
    pars.LRT.scaleksp=0;
    pars.LRT.niter = 5;
    LRT_nav_mask_enable = strcmp( pars.method, 'LRT');
    lrt_nav_weight = 2e-2;
    
    %--------------------------------------------------------------
    %% UPDATE mask
    
    
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
        sh_dim = length(recon_shot_range);
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
    subplot(211); montage(ttt, 'displayrange',[]);xlabel('kz'); ylabel('ky'); title('sampling mask shot-by-shot');
    subplot(223); imagesc(sum(ttt,4), [0 2]); colormap jet; axis off; axis equal; xlabel('kz'); ylabel('ky'); title('sample times'); colorbar;
    ttt = sum(bsxfun(@times, ttt, permute(1:size(ttt,4), [4 3 1 2])),4);
    subplot(224); imagesc(ttt, [0 sh_dim]); colormap jet; axis off; axis equal; xlabel('kz'); ylabel('ky'); title('sample shot #'); colorbar;
    clear ttt;
    
    %================separate navigator for LRT============================
    LRT_nav_mask = [];
    if(LRT_nav_mask_enable)
        LRT_nav_mask = zeros(size(mask));
        LRT_nav_mask(:,k0_idx(2)-10:k0_idx(2)+10,k0_idx(3)-8:k0_idx(3)+8,:,:) = 1;
        %         LRT_nav_mask = ones(size(mask));
    end
    
    %% SELECT Phase error map
    %select pe x loc
    
    %opt 1
    phase_error_3D_current =  phase_error_3D(recon_x_loc,:,:,:,:);
    
    %opt 2
    if(1)
        phase_error_3D_current =  permute(spiral_nav_im(recon_x_loc,:,:,:,:),[1 2 3 5 4]);
        phase_error_3D_current = bsxfun(@rdivide, phase_error_3D_current, phase_error_3D_current(:,:,:,:,1));
        phase_error_3D_current(find(isnan(phase_error_3D_current)+isinf(phase_error_3D_current))) = 0;
    end
    
    %opt3
    %     phase_error_3D_current = ones(size( phase_error_3D(recon_x_loc,:,:,:,:)));
    
    
    
    %% corrupt data, generate self-nav mask(POCS) and more
    
    %select channel ids
    
    if(length(ch_id)>1)
        sense_map_3D = normalize_sense_map(TSE_sense_map(recon_x_loc,:,:,ch_id));
    else
        sense_map_3D = ones(size(im_ref));
    end
    
    im_ref_ch =  bsxfun(@times, im_ref(recon_x_loc,:,:), sense_map_3D);
    im_b0_ref_ch =  bsxfun(@times, im_b0_ref(recon_x_loc,:,:), sense_map_3D);
    kspa_b0ref_ch = fft2d(squeeze(im_b0_ref_ch)); %b0 image has higher intensity
    s = ceil(TSE.kz_dim/2);
    %     figure(121); montage(permute(squeeze(abs(im_ref_ch)),[1 2 4 3]),'displayrange',[]); title('reference image, channel by channel');
    im_pe = bsxfun(@times,im_ref_ch,phase_error_3D_current);
    %======================================================================
    
    %================self_navigator size for POCS==========================
    if(strcmp( pars.method, 'POCS_ICE')) % for LRT separate navigator
        
        mask(:,k0_idx(2)-3:k0_idx(2)+3,k0_idx(3)-3:k0_idx(3)+3,:,:) = 1;
    end
    
    
    %================SENSE undersampling===================================
    if(simulate_ky_sense_2)
        mask(:,1:2:end,:,:,:) = 0;
    end
    if(simulate_kz_sense_2)
        mask(:,:,1:2:end,:,:) = 0;
    end
    
    
    kspa_xyz = bsxfun(@times, fft3d(im_pe), mask);
    
    
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
    
    
    %direct recon
    k_combine_shots=sum(kspa_xyz,5)./sum(abs(kspa_xyz)>0,5); %non zeros average
    k_combine_shots(find(isnan(k_combine_shots)))=0; k_combine_shots(find(isinf(k_combine_shots)))=0;
    
    im_ch_by_ch=ifft3d(k_combine_shots);
    im_recon_direct(recon_x_loc, :,:)=sqrt(sum(abs(im_ch_by_ch).^2, 4));
    %     figure(4); montage(permute(abs(im_recon_direct),[1 2 4 3]),'displayrange',[]); title('direct recon');
    clear  im_ch_by_ch k_combine_shots
    
    
    
    %% LRT navigator simulation
    spiral_nav = 1;
    
    %---- for LRT separate navigators --------
    if(~isempty(LRT_nav_mask))
        
        
        if(spiral_nav)  %corrupt the navigator or load a corrupted nav
            disp('>>>>>>>>>>>>>>>>>>>Spiral Navigator is generated based on images with phase error (im_pe) and trajectory!<<<<<<<<<<<<<<<<<<<<<<<')
            trj_scale = 2*pi/max((max(trj_simulation_kx)-min(trj_simulation_kx)),(max(trj_simulation_ky)-min(trj_simulation_ky)));
            trj_nufft = [trj_simulation_kx trj_simulation_ky] * trj_scale;
            [~, ky, kz, nc, nshot] = size(im_pe);
            A_nav = nuFTOperator(trj_nufft,[ky kz],ones(ky, kz),6);
            for c = 1:nc
                kspa_sosnav_b0ref_ch(:,c) = A_nav * squeeze(im_b0_ref_ch(1,:,:,c)); 
                for s = 1:nshot
                    kspa_x_yz_sosnav(:,c, s)= A_nav * squeeze(im_pe(1,:,:,c,s));
                end
            end
            A_nav_sense = nuFTOperator(trj_nufft,[ky kz],squeeze(sense_map_3D),6);
            disp('recon spiral nav....');
            for s=1:nshot
                %nav_im(:,:,1,s) = A_nav_sense' * col(kspa_x_yz_sosnav(:,:,s));
                nav_im(:,:,1,s) = regularizedReconstruction(A_nav_sense, col(kspa_x_yz_sosnav(:,:,s)),@L2Norm,0.1,'maxit',5);
            end
            disp('Done');
        else
            disp('>>>>>>>>>>>>>>>>>>>Cartesian Navigator is generated based on images with phase error (im_pe) and LRT_nav_mask!<<<<<<<<<<<<<<<<<<<<<<<')
            kspa_full = fft3d(im_pe);
            kspa_xyz_nav =  bsxfun(@times, kspa_full, LRT_nav_mask);
            
            kspa_x_yz_nav = ifft1d(kspa_xyz_nav);
            nav_im = ifft2d(squeeze(kspa_x_yz_nav));
            
            %calc b0 nav and kspace as well
            kspa_nav_b0ref_ch = kspa_b0ref_ch.*squeeze(abs(kspa_x_yz_nav(:,:,:,:,1))>0);
        end
            
        %----disp----%
        figure(411);
        subplot(131);montage(permute(squeeze(abs(nav_im(:,:,1,:))),[1 2 4 3]),'displayrange',[]); title('all nav images mag.')
        subplot(132);imshow(abs(squeeze(im_pe(1,:,:,1,1))),[]);
        subplot(133);montage(permute(squeeze(abs(LRT_nav_mask)),[1 2 4 3]),'displayrange',[]); title('mask')
        figure(412);
        nav_im_diff = bsxfun(@rdivide, nav_im, nav_im(:,:,:,1));
        nav_im_diff(find(isnan(nav_im_diff)+isinf(nav_im_diff)))=0;
        subplot(121); montage(permute(squeeze(angle(nav_im_diff(:,:,1,:))),[1 2 4 3]),'displayrange',[-pi pi]); colormap jet;  title('all nav images phase.')
        phase_error_3D_current_diff = bsxfun(@rdivide, phase_error_3D_current, phase_error_3D_current(:,:,:,:,1));
        phase_error_3D_current_diff(find(isnan(phase_error_3D_current_diff)+isinf(phase_error_3D_current_diff)))=0;
        subplot(122); montage(permute(squeeze(angle(phase_error_3D_current_diff)),[1 2 4 3]),'displayrange',[-pi pi]); colormap jet;  title('true phase error')
        
    end
    
    
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
    
    
    %=========select data. =====================================================
    sense_map = permute(sense_map_3D(1,:,:,:),[2 3 4 1]);
    phase_error = permute(squeeze(phase_error_3D_current(1,:,:,:, recon_shot_range)),[1 2 4 3]);
    kspa = permute(kspa_x_yz(1, :, :, :, recon_shot_range),[2 3 4 5 1]);
    %========================================================================
   
    %=====kspa for LRT is different=============================================
    clear kspa
    kspa_temp = permute(kspa_x_yz(1, :, :, :, recon_shot_range),[2 3 4 5 1]);
    %bart cs recon ref
    kspa_bart =  sum(kspa_temp,4)./sum(abs(kspa_temp)>0,4); kspa_bart(find(isnan(kspa_bart)))=0; kspa_bart(find(isinf(kspa_bart)))=0; 
    
    % for LRT, cat navigator
    if(strcmp( pars.method, 'LRT'))        
        if(spiral_nav) %spiral nav + cartesian images
            %%
            pars.LRT.mix_trajectory = 1;
            pars.LRT.NUFFT_nav_sense = nuFTOperator(trj_nufft,[size(sense_map,1) size(sense_map,2)],squeeze(sense_map),6);
            pars.LRT.NUFFT_nav_1ch = nuFTOperator(trj_nufft,[size(sense_map,1) size(sense_map,2)],ones(size(sense_map,1), size(sense_map,2)),6);
            pars.LRT.trj_length = length(trj_nufft);
            
            for idx1 =1:length(recon_shot_range)
                for idx2=1:2
                    if(idx2==1) %nav colume
                        kspa{idx1,idx2} =  kspa_x_yz_sosnav(:,:,recon_shot_range(idx2));
                    else %tse colume
                        kspa{idx1,idx2} =  kspa_temp(:,:,:,recon_shot_range(idx2));
                    end
                end
            end
            kspa_b0{1,1} =  kspa_sosnav_b0ref_ch;
            kspa_b0{1,2} =  kspa_b0ref_ch;
            kspa = cat(1, kspa_b0, kspa);
        else  %catesian nav
            kspa = cat(5, permute(kspa_x_yz_nav(1, :, :, :,  recon_shot_range) * lrt_nav_weight,[2 3 4 5 1]), kspa_temp);
            %cat "b0" info into the first row
            kspa = cat(4, cat(5, kspa_nav_b0ref_ch * lrt_nav_weight, kspa_b0ref_ch), kspa);
        end
    end
    
    
    
    
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
% figure(109);
% % subplot(141);
% montage(permute(abs(im_ref(recon_x_range,:,:)),[1 2 4 3]),'displayrange',[]); title('reference');
% figure(110);
% %     subplot(142);
% montage(permute(abs(im_recon_direct(recon_x_range,:,:)),[1 2 4 3]),'displayrange',[]); title('direct recon');
% figure(111);
% %     subplot(143);
% montage(permute(abs(image_corrected(recon_x_range,:,:)),[1 2 4 3]),'displayrange',[]); title('msDWIrecon');

toc;