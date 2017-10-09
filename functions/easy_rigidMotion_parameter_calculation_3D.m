function easy_rigidMotion_parameter_calculation_3D(fn)

    load(fn);


    ref_shot_idx = 4;
    [kx, ky, kz, ch, shot] = size(nav_kspa_recon_nufft);

    nav_im_recon_nufft_diff = bsxfun(@rdivide, nav_im_recon_nufft, nav_im_recon_nufft(:,:,:,:,1,:));
    nav_im_recon_nufft_diff(find(isinf(nav_im_recon_nufft_diff)))=0;
    nav_im_recon_nufft_diff(find(isnan(nav_im_recon_nufft_diff)))=0;
        
    for ch_idx = 1:ch

          ref_nav = nav_kspa_recon_nufft(:,:,:,ch_idx,ref_shot_idx); %ref_shot_idx shot as the reference navigator

  
        % --------- Linear &GLOBAL phase errror estimation in image domain
        


        min_kx = -ceil((kx-1)/2);
        min_ky = -ceil((ky-1)/2);
        min_kz = -ceil((kz-1)/2);

        [X, Y, Z] = meshgrid(min_kx:kx-1+min_kx, min_ky:ky-1+min_ky, min_kz:kz-1+min_kz); 
        fit_x(:,1) = X(:);
        fit_x(:,2) = Y(:);
        fit_x(:,3) = Z(:);
        for shot_idx = 1:shot

%              fitting_nav_ima = exp(i * (nav_ima_global_calibrated_phase_diff(:,:,ch_idx,shot_idx)));  %use cpx --> don't need to worry about phase wrapping

        %     %----------lsqcurvefit ------------- NOT WORKING WELL
        %     f = fitting_nav_ima(:);
        %     
        %     fun = @(c, x)exp(i * (fit_x(:,1) * c(1) + fit_x(:,2) * c(2)));
        %     
        %     % solve with lsqcurvefit
        %     options = optimset('TolX',1e-6);
        %     c0 = [0 0];
        %     cc = lsqcurvefit(fun, c0, fit_x, f, [], [], options);
        %     lfit = fun(cc, fit_x);
        %     fitting_result(:,:,ch_idx,shot_idx) = reshape(lfit, [kx ky]);

        %   %----------2D lsqnonlin------------- ...

            clear  fitting_nav_ima_phase_unwrapped idx_id fun_b options b0 b fitting_results_1d fitting_results_2d
            fitting_nav_ima_raw_phase_diff = angle(squeeze(nav_im_recon_nufft_diff(:,:,:,ch_idx,shot_idx)));
 
            all_intensity = sort(col(abs(nav_im_recon_nufft(:,:,:,ch_idx,shot_idx))));
            
            % pixels with high magnitude; not NaN; not Inf; being used for fitting
            image_intensity_cutoff = all_intensity(round(length(all_intensity)*0.5)); %0.5: 50% of the pixels were used; 0.8: 20% of the pixels were used

            idx_id = find((abs(nav_im_recon_nufft(:,:,:,ch_idx,shot_idx))>image_intensity_cutoff)&(~isnan(fitting_nav_ima_raw_phase_diff))&(~isinf(fitting_nav_ima_raw_phase_diff)));
            [sub_x, sub_y, sub_z] = ind2sub([kx  ky  kz], idx_id );
            
            fun_b = @(b)abs(   exp(1i*(   sub_x * b(1) + sub_y * b(2) +  sub_z * b(3) +b(4) ) ) ...
                          -  exp(1i*fitting_nav_ima_raw_phase_diff(idx_id)));

            options.Algorithm = 'trust-region-reflective';
            options.MaxFunEvals = 1000;
            b0 = [0 0 0 0];
            b = lsqnonlin(fun_b,b0,[], [], options);

            fitting_results_1d = zeros(kx * ky * kz, 1);
            fitting_results_1d(idx_id) = sub_x * b(1) + sub_y * b(2) +  sub_z * b(3) +b(4); 
            fitting_results_3d = reshape(fitting_results_1d, kx, ky, kz);
            figure(11); 
            subplot(121);immontage4D(fitting_nav_ima_raw_phase_diff,[ -5*pi 5*pi]); colormap jet
            title(['shot: ', num2str(shot_idx)]);
            subplot(122);immontage4D(fitting_results_3d,[ -5*pi 5*pi]); colormap jet
            title(['ch: ', num2str(ch_idx)]);
            drawnow();
            fitted_nav_ima_phase(:,:,:,ch_idx,shot_idx) = fitting_results_3d;
            linear_phase_xyz(:,ch_idx,shot_idx) = b(1:3); %in rad per navigator pixel
            global_phase(ch_idx,shot_idx) = b(4); %in rad

        %   %----------1D lsqnonlin------------- SEEMS GOOD; but how to combine to dirrections???
            %====================X
        %     fitting_nav_ima_1dx = fitting_nav_ima(:,12);
        % 
        %     idx_1d = find(abs(nav_ima_global_calibrated(:,12,ch_idx,1))>1e6);
        %     idx_1d = idx_1d';
        %     i = sqrt(-1);
        %     % with wrapping issue
        % %             fun_c_x = @(c_x)(abs(  exp(i * (idx_1d * c_x(1) + c_x(2)))  -  fitting_nav_ima_1dx(idx_1d)  ).^2);
        %     % consider wrapping issue    
        % %             fun_c_x = @(c_x)(  (idx_1d * c_x(1) + c_x(2) )  -  unwrap(angle(fitting_nav_ima_1dx(idx_1d)))  );
        % 
        %     options.Algorithm = 'trust-region-reflective';
        %     options.MaxFunEvals = 1000;
        % %     options = optimset('TolX',1e-8);
        %     c0 = [0 global_phase_diff_initial(shot_idx)];
        %     c_x = lsqnonlin(fun_c_x,c0,[], [], options)
        % 
        %     fitting_results = idx_1d * c_x(1) + c_x(2);
        %      
        %     %====================Y
        %     fitting_nav_ima_1dy = fitting_nav_ima(12,:);
        %     
        %     idy_1d = find(abs(nav_ima_global_calibrated(12,:,ch_idx,1))>1e6);
        %     i = sqrt(-1);
        %     % with wrapping issue
        % %     fun_c_x = @(c_x)abs(  exp(i * (idx_1d * c_x ))  -  fitting_nav_ima_1dx(idx_1d)  );
        %     % consider wrapping issue    
        %     fun_c_y = @(c_y)(  (idy_1d * c_y(1) + c_y(2) )  -  unwrap(angle(fitting_nav_ima_1dy(idy_1d)))  );
        %     
        %     options.Algorithm = 'trust-region-reflective';
        %     options.MaxFunEvals = 1000;
        % %     options = optimset('TolX',1e-8);
        %     c0 = [10 10];
        %     c_y = lsqnonlin(fun_c_y,c0,[], [], options)
        %     
        %     fitting_results = idy_1d * c_y(1) + c_y(2);
        %   
        %             

        end

        % -----DISPLAY

          figure(7); montage(nav_ima_phase_unwrapped_diff(:,:,ch_idx,:),'displayrange',[-5*pi 5*pi]); colormap jet
          figure(8); montage(fitted_nav_ima_phase(:,:,ch_idx,:),'displayrange',[-5*pi 5*pi]); colormap jet

    end
    save(fn, 'nav_ima_phase_unwrapped_diff','fitted_nav_ima_phase', 'linear_phase_xy','global_phase','global_phase_diff_initial','-append');
end