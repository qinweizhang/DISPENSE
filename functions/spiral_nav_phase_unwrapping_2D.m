function nav_ima_phase_unwrapped = spiral_nav_phase_unwrapping_2D(nav_im_to_unwrap, ref_shot_ix)
%INPUT
%
%nav_im_to_unwrap: complex data in size of [nav_dimx, nav_dimy, nav_chan, nav_shot] 
%
%OUTPUT
%
%nav_ima_phase_unwrapped: phase images in size of [nav_dimx, nav_dimy, nav_chan, nav_shot] 


[nav_dimx, nav_dimy, nav_chan, nav_shot] = size(nav_im_to_unwrap)



%------------------------------------------------------PHSE UNWRAPPING----------------------------------%

Manually_seed_point = 0;
for ch_idx = 1:nav_chan
    for shot_idx = 1:nav_shot
        nav_ima_phase_unwrapped(:,:,ch_idx,shot_idx) = GoldsteinUnwrap2D_r1_func(nav_im_to_unwrap(:,:,ch_idx,shot_idx), Manually_seed_point);
        figure(22);
        subplot(121); imshow(squeeze(angle(nav_im_to_unwrap(:,:,ch_idx,shot_idx))),[-pi pi]); title('L: wraped R: unwrapped ') ;  colorbar
        subplot(122); imshow(squeeze(nav_ima_phase_unwrapped(:,:,ch_idx,shot_idx)),[-pi pi]); colormap jet;title(['shot: ',num2str(shot_idx),'ch: ',num2str(ch_idx)]);  colorbar
    end
end



ch_idx = 1;
% ===================display wrapped phase
phase_wrapped = angle(nav_im_to_unwrap);

figure(601); title('display wrapped phase');

montage((phase_wrapped(:,:,ch_idx,:)),'displayrange',[-pi pi]); colormap jet; title(['all shots from ch', num2str(ch_idx)]);
% ===================display wrapped phase difference

phase_test_ch = phase_wrapped(:,:,ch_idx,:);
for shot_idx = 1:nav_shot
    phase_test_ch_diff(:,:,1,shot_idx) = phase_test_ch(:,:,1,shot_idx) -phase_test_ch(:,:,1,ref_shot_ix) ;
end
figure(602); montage(phase_test_ch_diff,'displayrange',[-pi pi]); colormap jet; title(['all shots phase differece from ch', num2str(ch_idx)]);


% ====================display unwrapped phase difference
phase_test_ch = nav_ima_phase_unwrapped(:,:,ch_idx,:);
for shot_idx = 1:nav_shot
    phase_test_ch_diff(:,:,1,shot_idx) = phase_test_ch(:,:,1,shot_idx) -phase_test_ch(:,:,1,ref_shot_ix) ;
end
figure(604); montage(phase_test_ch_diff,'displayrange',[-5*pi 5*pi]); colormap jet; title(['unwrapped phase difference:all shots phase differece from ch', num2str(ch_idx)]);

end