function[ nav_ima_phase_unwrapped_diff,fitted_nav_ima_phase, linear_phase_xy,global_phase] =  easy_rigidMotion_parameter_calculation( nav_ima_phase_unwrapped, nav_im_recon_nufft)

    ref_shot_idx = 4;
    [kx, ky, ch, shot] = size(nav_ima_phase_unwrapped);

    for ch_idx = 1:ch
        

%         figure(3); montage(abs(nav_kspa_recon_nufft(:,:,ch_idx,:)),'displayrange',[])
%         figure(4); montage(angle(nav_kspa_recon_nufft(:,:,ch_idx,:)),'displayrange',[-pi pi]); colormap jet
%         ref_nav = nav_kspa_recon_nufft(:,:,ch_idx,ref_shot_idx); %ref_shot_idx shot as the reference navigator
        
        if(0)
            for shot_idx = 1:shot
                
                current_nav = nav_kspa_recon_nufft(:,:,ch_idx,shot_idx);
                
                % argmax|Sref(K)|
                Ref_max_ix = find(abs(ref_nav(:)) == max(abs(ref_nav(:))));
                Ref_max_kx = (mod(mod(Ref_max_ix,(kx*ky)),kx));
                Ref_max_ky = (ceil(mod(Ref_max_ix,(kx*ky))/kx));
                Ref_max = ref_nav(Ref_max_kx, Ref_max_ky);
                % argmax|Sb(K)|
                current_max_ix = find(abs(current_nav(:)) == max(abs(current_nav(:))));
                current_max_kx = (mod(mod(current_max_ix,(kx*ky)),kx));
                current_max_ky = (ceil(mod(current_max_ix,(kx*ky))/kx));
                current_max = current_nav(current_max_kx, current_max_ky);
                
                global_phase_diff_initial(shot_idx) = angle(current_max) - angle(Ref_max); %rad
                
                nav_kspa_global_calibrated(:,:,ch_idx,shot_idx) = current_nav.*exp(-i * global_phase_diff_initial(shot_idx));
                
                offset_x(shot_idx) = current_max_kx - Ref_max_kx; %in pixel
                offset_y(shot_idx) = current_max_ky - Ref_max_ky; %in pixel
                
            end
            
            
            
            %------------------------------------------------------DISPLAY----------------------------------%
            
            
            
            for shot_idx = 1:shot
                nav_ima_raw_phase_diff(:,:,ch_idx,shot_idx) = angle(nav_im_recon_nufft(:,:,ch_idx,shot_idx)) - angle(nav_im_recon_nufft(:,:,ch_idx,ref_shot_idx));
            end
            figure(5); montage(nav_ima_raw_phase_diff(:,:,ch_idx,:),'displayrange',[-pi pi]); colormap jet
            
            for shot_idx = 1:shot
                nav_ima_phase_unwrapped_diff(:,:,ch_idx,shot_idx) = nav_ima_phase_unwrapped(:,:,ch_idx,shot_idx) - nav_ima_phase_unwrapped(:,:,ch_idx,ref_shot_idx);
            end
            figure(6); montage(nav_ima_phase_unwrapped_diff(:,:,ch_idx,:),'displayrange',[-pi pi]); colormap jet
            
            nav_ima_global_calibrated = bart('fft -i 3', nav_kspa_global_calibrated);
            for shot_idx = 1:shot
                nav_ima_global_calibrated_phase_diff(:,:,ch_idx,shot_idx) = angle(nav_ima_global_calibrated(:,:,ch_idx,shot_idx)) - angle(nav_ima_global_calibrated(:,:,ch_idx,ref_shot_idx));
            end
            figure(7); montage(nav_ima_global_calibrated_phase_diff(:,:,ch_idx,:),'displayrange',[-pi pi]); colormap jet
            
            %--------------------------------------------------------------------------------------------------%
            
            % --------- Linear &GLOBAL phase errror estimation in image domain
        end
        
        
        for shot_idx = 1:shot
            nav_ima_phase_unwrapped_diff(:,:,ch_idx,shot_idx) = nav_ima_phase_unwrapped(:,:,ch_idx,shot_idx) - nav_ima_phase_unwrapped(:,:,ch_idx,ref_shot_idx);
        end
        figure(6); montage(nav_ima_phase_unwrapped_diff(:,:,ch_idx,:),'displayrange',[-pi pi]); colormap jet
        
        
        [X, Y] = meshgrid(-(kx/2-1):kx/2, -(ky/2-1):ky/2);
        fit_x(:,1) = X(:);
        fit_x(:,2) = Y(:);
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
            fitting_nav_ima_phase_unwrapped = squeeze(nav_ima_phase_unwrapped_diff(:,:,ch_idx,shot_idx));
 
            all_intensity = sort(col(abs(nav_im_recon_nufft(:,:,ch_idx,shot_idx))));
            
            % pixels with high magnitude; not NaN; not Inf; being used for fitting
            image_intensity_cutoff = all_intensity(round(length(all_intensity)*0.5)); %0.5: 50% of the pixels were used; 0.8: 20% of the pixels were used

            idx_id = find((abs(nav_im_recon_nufft(:,:,ch_idx,shot_idx))>image_intensity_cutoff)&(~isnan(fitting_nav_ima_phase_unwrapped))&(~isinf(fitting_nav_ima_phase_unwrapped)));
            fun_b = @(b)(   (   (mod(idx_id-1,kx)-kx/2) * b(1) + (floor((idx_id-1)/kx)-ky/2) * b(2) +  b(3)  ) ...
                          -  fitting_nav_ima_phase_unwrapped(idx_id));

            options.Algorithm = 'trust-region-reflective';
            options.MaxFunEvals = 1000;
%             b0 = [0 0 global_phase_diff_initial(shot_idx)];
            b0 = [0 0 0];

            b = lsqnonlin(fun_b,b0,[], [], options);

            fitting_results_1d = zeros(kx * ky, 1);
            fitting_results_1d(idx_id) =  (mod(idx_id-1,kx)-(kx/2)) * b(1) + (floor((idx_id-1)/kx)-(ky/2)) * b(2) +  b(3); 
            fitting_results_2d = reshape(fitting_results_1d, kx, ky);
            figure(11); 
            subplot(121);imshow(fitting_nav_ima_phase_unwrapped,[ -5*pi 5*pi]); colormap jet
            title(['shot: ', num2str(shot_idx)]);
            subplot(122);imshow(fitting_results_2d,[ -5*pi 5*pi]); colormap jet
            title(['ch: ', num2str(ch_idx)]);
            drawnow();
            fitted_nav_ima_phase(:,:,ch_idx,shot_idx) = fitting_results_2d;
            linear_phase_xy(:,ch_idx,shot_idx) = b(1:2); %in rad per navigator pixel
            global_phase(ch_idx,shot_idx) = b(3); %in rad

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
%     save(fn, 'nav_ima_phase_unwrapped_diff','fitted_nav_ima_phase', 'linear_phase_xy','global_phase','global_phase_diff_initial','-append');
%         save(fn, 'nav_ima_phase_unwrapped_diff','fitted_nav_ima_phase', 'linear_phase_xy','global_phase','-append');

end