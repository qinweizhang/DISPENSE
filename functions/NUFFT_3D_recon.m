% OUTPUT
%
% nav_im_recon_nufft: reconed navigator image in size of [x,y,z,ch,shot]
%
% 
% INPUT
% 
% nav_k_spa_data
% trj_fn
% recon_par
% nav_sense_map
% varargin{1} 
% varargin{2:3} 


function nav_im_recon_nufft = NUFFT_3D_recon(nav_k_spa_data,trj_fn,recon_par,nav_sense_map, varargin)


load(trj_fn);
if( nargin==5 )
    nav_sense_Psi = double(varargin{1});
elseif(nargin >5 )
    assert(nargin == 7, 'Input both Offcenter and FOV in mm');
    nav_sense_Psi = varargin{1};
    Offcenter_xy = varargin{2};
    FOV_xy = varargin{3};
end

%% display singal quality

dist = abs(trj_meas_kx + i .* trj_meas_ky);
dist = abs(dist + i.*trj_meas_kz);
figure(16); plot(dist./max(abs(dist(:))));
hold on
plot((abs(nav_k_spa_data(:,2,1,1))./max(abs(nav_k_spa_data(:,2,1,1)))),'r')
title('signal vs time points');legend('trajectory distance','abs(signal)')

selected_point = 1:length(trj_meas_kx);
if(isfield(recon_par, 'end_point') && isfield(recon_par, 'skip_point'))
    if(isempty(recon_par.end_point))
        recon_par.end_point = length(trj_meas_kx);
    end 
    selected_point = recon_par.skip_point+1:recon_par.end_point;
end
if(isfield(recon_par, 'selected_point') )
     if(isempty(recon_par.selected_point))
        recon_par.selected_point = 1:length(trj_meas_kx);
    end 
    selected_point = recon_par.selected_point;
end
%% Scale trajectory

trj_meas_kx_t = squeeze(trj_meas_kx(selected_point,1));
trj_meas_ky_t = squeeze(trj_meas_ky(selected_point,1));
trj_meas_kz_t = squeeze(trj_meas_kz(selected_point,1));

clear trj_meas_kx trj_meas_ky trj_meas_kz k_spa_data
trj_meas_kx = trj_meas_kx_t';
trj_meas_ky = trj_meas_ky_t';
trj_meas_kz = trj_meas_kz_t';

clear trj_meas_kx_t trj_meas_ky_t trj_meas_kz_t sig_bart_t

try
    recon_scale_factor = recon_par.recon_dim./ recon_par.acq_dim;
catch
    recon_scale_factor = [1 1 1];
end

scale_foctor_xy = max(2*pi/(max(trj_meas_kx)-min(trj_meas_kx))./recon_scale_factor(1),2*pi/(max(trj_meas_ky)-min(trj_meas_ky))./recon_scale_factor(2));
trj_meas_kx_scaled = trj_meas_kx * scale_foctor_xy;
trj_meas_ky_scaled = trj_meas_ky * scale_foctor_xy;

scale_foctor_z = 2*pi/(max(trj_meas_kz)-min(trj_meas_kz))./recon_scale_factor(3);

if(recon_par.ignore_kz == 1)
    scale_foctor_z = 0;
end
trj_meas_kz_scaled = trj_meas_kz * scale_foctor_z;


figure(8);
plot(trj_meas_kx_scaled); hold on; plot(trj_meas_ky_scaled,'r');  plot(trj_meas_kz_scaled,'k');legend('measured kx (a.u.)','measured ky (a.u.)','measured kz (a.u.)');

figure(801);
plot3(trj_meas_kx_scaled, trj_meas_ky_scaled, trj_meas_kz_scaled); title('-pi to pi')


trj_nufft = double(cat(1, trj_meas_kx_scaled, trj_meas_ky_scaled, trj_meas_kz_scaled));
trj_nufft = trj_nufft';
clear im_recon_nufft nav_im_recon_nufft

%% recon

% ========Orthogonal SENSE maps: recombine coils to make the Psi map indentity mtx (SNR optimized)
disp('Orthogonal SENSE maps and nav kspace data...');
tic
if(exist('nav_sense_Psi','var'))
    if(~isempty(nav_sense_Psi))
        for t=1:length(nav_sense_Psi)  %make sure diag(nav_sense_Psi) are all real; they should be!
            nav_sense_Psi(t,t) = real(nav_sense_Psi(t,t));
        end
        L = chol(nav_sense_Psi,'lower'); %Cholesky decomposition; L: lower triangle
        L_inv = inv(L);
        for c = 1:size(nav_sense_Psi,1)
            %recombine sense map
            sense_map_orthocoil(:,:,:,c) = sum(bsxfun(@times, nav_sense_map, permute((L_inv(c,:)),[1 3 4 2])), 4);
            
            %recombine kspa
            nav_k_spa_data_orthocoil(:,c,:,:) = sum(bsxfun(@times, nav_k_spa_data(:,:,:,:), (L_inv(c,:))), 2);
        end
%         sense_map_orthocoil = nav_sense_map;  %rewrite 
        
        figure(401);
        s = round(size(nav_sense_map,3)/2);
        subplot(121); montage(abs(nav_sense_map(:,:,s,:)),'displayrange',[]); title('originial SENSE map')
        subplot(122); montage(abs(sense_map_orthocoil(:,:,s,:)),'displayrange',[]); title('Orthogonal SENSE map')
        
        nav_sense_map =  sense_map_orthocoil;
        clear kspa_orthocoil nav_k_spa_data_orthocoil
        %renormalize sense
        nav_sense_map = normalize_sense_map(nav_sense_map);
    end
end
toc
% =================================================================================================



if(recon_par.recon_all_shot)
    end_shot_idx = size(nav_k_spa_data,3);
else
    end_shot_idx = 1;
end

if(recon_par.parfor)  %parfor recon
    parfor shot_nr = 1: end_shot_idx
        shot_nr
        
        sig_kspa = double(nav_k_spa_data(selected_point,:,shot_nr,recon_par.dyn_nr));
        
        %============xy offset compensation: based on trajectory======================================================
        %if(exist('Offcenter_xy','var'))  %error in parfor
            %---x offset
            if(Offcenter_xy(1)~=0)
                ima_offcenter_FOV_ratio = Offcenter_xy(1)/FOV_xy(1);
                kspa_linear_phase_rate = ima_offcenter_FOV_ratio * (2*pi);                  % in (rad/kspce pixel)
                kspa_trj_in_pixel = trj_nufft(:,1) / pi * (recon_par.recon_dim(1) * 0.5);   %convert trajectory unit from rad to pixel
                kspa_clibration_phase = kspa_trj_in_pixel * kspa_linear_phase_rate;
                sig_kspa = bsxfun(@times, sig_kspa, exp(i*kspa_clibration_phase));          %add this calibration linear phase
            end
            %---y offset
            if(Offcenter_xy(2)~=0)
                ima_offcenter_FOV_ratio = Offcenter_xy(2)/FOV_xy(1);                  %still use FOV_xy(1) as spiral FOV is always squared
                kspa_linear_phase_rate = ima_offcenter_FOV_ratio * (2*pi);                  % in (rad/kspce pixel)
                kspa_trj_in_pixel = trj_nufft(:,2) / pi * (recon_par.recon_dim(2) * 0.5);   %convert trajectory unit from rad to pixel
                kspa_clibration_phase = kspa_trj_in_pixel * kspa_linear_phase_rate;
                sig_kspa = bsxfun(@times, sig_kspa, exp(i*kspa_clibration_phase));          %add this calibration linear phase
                
            end
        %end
        %===========================================================================================================
        
        
        if(~recon_par.channel_by_channel) %(recon_par.sense_map_recon) %all in one recon
                        
            sig_nufft = col(double(sig_kspa));
            
            A=nuFTOperator(trj_nufft,recon_par.recon_dim,double(nav_sense_map),6);
            
            % simple inverse
            nav_im_recon_nufft_direct_inverse = A'*sig_nufft;
            figure(702); montage(permute(abs(nav_im_recon_nufft_direct_inverse), [1 2 4 3]),'displayrange',[]); title('direct inverse');
            
            
            %call CG-SENSE with L2-norm regularization
            nav_im_recon_nufft(:,:,:,:,shot_nr)=regularizedReconstruction(A,sig_nufft,@L2Norm,recon_par.lamda,'maxit',recon_par.interations,'tol', 1e-10);
            %                 nav_im_recon_nufft(:,:,:,:,shot_nr)=regularizedReconstruction(A,sig_nufft,'maxit',recon_par.interations);
            
        end
    end
else %no parfor recon
    for shot_nr = 1: end_shot_idx
        shot_nr
        
        sig_kspa = double(nav_k_spa_data(selected_point,:,shot_nr,recon_par.dyn_nr));
        
        %============xy offset compensation: based on trajectory======================================================
        if(exist('Offcenter_xy','var'))
            %---x offset
            if(Offcenter_xy(1)~=0)
                ima_offcenter_FOV_ratio = Offcenter_xy(1)/FOV_xy(1);
                kspa_linear_phase_rate = ima_offcenter_FOV_ratio * (2*pi);                  % in (rad/kspce pixel)
                kspa_trj_in_pixel = trj_nufft(:,1) / pi * (recon_par.recon_dim(1) * 0.5);   %convert trajectory unit from rad to pixel
                kspa_clibration_phase = kspa_trj_in_pixel * kspa_linear_phase_rate;
                sig_kspa = bsxfun(@times, sig_kspa, exp(i*kspa_clibration_phase));          %add this calibration linear phase
            end
            %---y offset
            if(Offcenter_xy(2)~=0)
                ima_offcenter_FOV_ratio = Offcenter_xy(2)/FOV_xy(1);                  %still use FOV_xy(1) as spiral FOV is always squared
                kspa_linear_phase_rate = ima_offcenter_FOV_ratio * (2*pi);                  % in (rad/kspce pixel)
                kspa_trj_in_pixel = trj_nufft(:,2) / pi * (recon_par.recon_dim(2) * 0.5);   %convert trajectory unit from rad to pixel
                kspa_clibration_phase = kspa_trj_in_pixel * kspa_linear_phase_rate;
                sig_kspa = bsxfun(@times, sig_kspa, exp(i*kspa_clibration_phase));          %add this calibration linear phase
                
            end
        end
        %===========================================================================================================
        
        
        if(~recon_par.channel_by_channel) %(recon_par.sense_map_recon) %all in one recon
            
            
            sig_nufft = col(double(sig_kspa));
            
            A=nuFTOperator(trj_nufft,recon_par.recon_dim,double(nav_sense_map),6);
            
            % simple inverse
            nav_im_recon_nufft_direct_inverse = A'*sig_nufft;
            figure(702); montage(permute(abs(nav_im_recon_nufft_direct_inverse), [1 2 4 3]),'displayrange',[]); title('direct inverse');
            
            
            %call CG-SENSE with L2-norm regularization
            nav_im_recon_nufft(:,:,:,:,shot_nr)=regularizedReconstruction(A,sig_nufft,@L2Norm,recon_par.lamda,'maxit',recon_par.interations,'tol', 1e-10);
            %                 nav_im_recon_nufft(:,:,:,:,shot_nr)=regularizedReconstruction(A,sig_nufft,'maxit',recon_par.interations);
            
        else %ch by ch recon
            
            disp('recon ch by ch!')
            
            for ch =1:size(nav_k_spa_data,2)
                sig_nufft = double(sig_kspa(:, ch));
                sig_nufft = sig_nufft';
                A=nuFTOperator(trj_nufft,recon_par.recon_dim,nav_sense_map(:,:,:,ch),6);
                
                % simple inverse
                %                 nav_im_recon_nufft(:,:,:,ch) = A'*sig_nufft';
                
                
                
                %call CG-SENSE with L2-norm regularization
                nav_im_recon_nufft(:,:,:,ch,shot_nr)=regularizedReconstruction(A,sig_nufft',@L2Norm,recon_par.lamda,'maxit',recon_par.interations);
                %             nav_im_recon_nufft(:,:,:,ch,shot_nr)=regularizedReconstruction(A,sig_nufft','maxit',recon_par.interations);
            end
            % im_recon_nufft = flipdim(flipdim(im_recon_nufft,1),2);
        end
    end
end


%% display
display_shot_nr = 1;
try
    if(recon_par.sense_map_recon) %all in one recon
        
        figure(37);
        subplot(121); montage(permute(squeeze(abs(nav_im_recon_nufft)), [1 2 4 3]),'displayrange',[]); title('all slices rss'); xlabel('P'); ylabel('M')
        subplot(122); montage(permute(squeeze(abs(nav_im_recon_nufft)), [2 3 4 1]),'displayrange',[]); title('all slices rss'); xlabel('S'); ylabel('P')
        
        figure(38);
        subplot(121); montage(permute(squeeze(angle(nav_im_recon_nufft)), [1 2 4 3]),'displayrange',[-pi pi]); title('phase all slices'); xlabel('P'); ylabel('M'); colormap jet;
        subplot(122); montage(permute(squeeze(angle(nav_im_recon_nufft)), [2 3 4 1]),'displayrange',[-pi pi]); title('phase all slices'); xlabel('S'); ylabel('P'); colormap jet;
        
    else
        
        figure(35); montage(permute(abs(squeeze(nav_im_recon_nufft(:,:,:,6,display_shot_nr))), [1 2 4 3]),'displayrange',[]); title('all slices ch6')
        figure(36); montage(-1* permute(angle(squeeze(nav_im_recon_nufft(:,:,:,6,display_shot_nr))), [1 2 4 3]),'displayrange',[-pi pi]); colormap jet; title('all slices ch6')
        
        nav_im_recon_nufft_rss = squeeze(bart('rss 8', nav_im_recon_nufft));
        figure(37);
        subplot(121); montage(permute(abs(nav_im_recon_nufft_rss), [1 2 4 3]),'displayrange',[]); title('all slices rss'); xlabel('P'); ylabel('M')
        subplot(122); montage(permute(abs(nav_im_recon_nufft_rss), [2 3 4 1]),'displayrange',[]); title('all slices rss'); xlabel('S'); ylabel('P')
        
        figure(38);
        subplot(121); immontage4D(permute(squeeze(angle(nav_im_recon_nufft)), [1 2 4 3]), [-pi pi]); title('phase all slices'); xlabel('P'); ylabel('M'); colormap jet; colorbar
        subplot(122); immontage4D(permute(squeeze(angle(nav_im_recon_nufft)), [2 3 4 1]), [-pi pi]); title('phase all slices'); xlabel('S'); ylabel('P'); colormap jet; colorbar
        
    end
catch
    warning('Kerry:  figure display not succeed!');
    
end
end