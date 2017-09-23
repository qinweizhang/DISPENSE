size(phase_error_3D)

phase_error_3D_yzx = permute(phase_error_3D, [2 3 1 4]);

figure(61); slice = 100;
montage(angle(phase_error_3D_yzx(:,:,slice,:)),'displayrange',[-pi pi]); colormap jet

tt = reshape(phase_error_3D_yzx(:,:,slice,:),size(phase_error_3D_yzx,1)*size(phase_error_3D_yzx,2),size(phase_error_3D_yzx,4));
[U S V] = svd(squeeze(tt));  figure(62); plot(abs(diag(S)))

diag_S = diag(S);
rank = 7;
S_lr_s = diag([diag_S(1:rank);zeros(length(diag_S)-rank,1)]);
S_lr = cat(1,S_lr_s,zeros(size(S,1)-size(S,2),size(S,2)));
tt_low_rank = U*S_lr*V';

phase_error_3D_yzx_lr = reshape(tt_low_rank, size(phase_error_3D_yzx,1), size(phase_error_3D_yzx,2),1, size(phase_error_3D_yzx,4) );

figure(63)
montage(angle(phase_error_3D_yzx_lr),'displayrange',[-pi pi]); colormap jet

%%
ref_shot =1; slice = 100;
phase_error_3D_yzx_unwrapped = spiral_nav_phase_unwrapping_2D(phase_error_3D_yzx(:,:,slice,10),ref_shot);
figure(62); 

subplot(121); montage(angle(phase_error_3D_yzx(:,:,slice,:)),'displayrange',[-pi pi]); colormap jet
subplot(122); montage(phase_error_3D_yzx_unwrapped,'displayrange',[-pi pi]); colormap jet