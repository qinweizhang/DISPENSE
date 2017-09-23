function nav_im_recon_nufft = NUFFT_3D_recon(nav_k_spa_data,trj_fn,recon_par,nav_sense_map)

%OUTPUT
%
% nav_im_recon_nufft: reconed navigator image in size of [x,y,z,ch,shot]

load(trj_fn); 

% display singal quality
dist = abs(trj_meas_kx + i .* trj_meas_ky);
dist = abs(dist + i.*trj_meas_kz);
figure(16); plot(dist./max(abs(dist(:))));
hold on
plot((abs(nav_k_spa_data(:,2,1,1))./max(abs(nav_k_spa_data(:,2,1,1)))),'r')
title('signal vs time points');legend('trajectory distance','abs(signal)')


if(isempty(recon_par.end_point))
    recon_par.end_point = length(trj_meas_kx);
end

selected_point = recon_par.skip_point+1:recon_par.end_point;




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
if(recon_par.recon_all_shot)
    end_shot_idx = size(nav_k_spa_data,3);
else
    end_shot_idx = 1;
end
for shot_nr = 1: end_shot_idx
    shot_nr
    
    sig_kspa = nav_k_spa_data(selected_point,:,shot_nr,recon_par.dyn_nr);
    
    if(~recon_par.channel_by_channel) %(recon_par.sense_map_recon) %all in one recon
        sig_nufft = col(double(sig_kspa));
        
%         nav_sens_map = conj(nav_sens_map);
        A=nuFTOperator(trj_nufft,recon_par.recon_dim,nav_sense_map,6);
        
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