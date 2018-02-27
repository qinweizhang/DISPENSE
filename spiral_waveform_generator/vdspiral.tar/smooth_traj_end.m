function k_optimized = smooth_traj_end(k)

gr_dwell    = 6.4e-6;   %in s
gamma       = 4257.6;

phi_end = angle(k(end));
k =  k .* exp(-i.*phi_end); %always ends at 0 degree

g_y = 1/gamma * cat(2, 0, diff(imag(k))./gr_dwell);

g_x = 1/gamma * cat(2, 0, diff(real(k))./gr_dwell);

fix_gr_end_str = false;
fix_gr_end_angle = true;
brutal = false;
assert(sum([fix_gr_end_str fix_gr_end_angle brutal])==1);

if(fix_gr_end_str)
    %% control the strength of the end gradient to (0)
    
    %put g_x(end) g_y(end) to res * g_end
    blip_sr = 9000; %G/cm/s  90 mT/m/ms
    res = 0;
    
    slop = (1-res)*max(abs(g_x(end)/blip_sr), abs(g_y(end)/blip_sr)); %seconds
    slop = round(slop/gr_dwell)*gr_dwell;
    blip_sr_x = -g_x(end)/slop;
    blip_sr_y = -g_y(end)/slop;
    
    str_step_x = blip_sr_x * gr_dwell;
    g_x_tail = [g_x(end):str_step_x:res*g_x(end)];
    g_x_tail = g_x_tail(2:end);
    
    str_step_y = blip_sr_y * gr_dwell;
    g_y_tail = [g_y(end):str_step_y:res*g_y(end)];
    g_y_tail = g_y_tail(2:end);
    
    g_x_final = [g_x g_x_tail];
    g_y_final = [g_y g_y_tail];
    
    
    kx = zeros(1, length(g_x_final));
    ky = zeros(1, length(g_x_final));
    for idx = 2:length(g_x_final)
        kx(idx) = kx(idx-1) + gamma*gr_dwell*g_x_final(idx);
        ky(idx) = ky(idx-1) + gamma*gr_dwell*g_y_final(idx);
    end
    k_optimized = kx + i*ky;
    
    phi_end = angle(k_optimized(end));
    k_optimized =  k_optimized .* exp(-i.*phi_end); %always ends at 0 degree
    disp(['Smoothing tailor length: ', num2str(length(g_y_tail))]);
elseif(fix_gr_end_angle)
    %% control the phase of the end gradient to be perpendicular to k(end): g(end)*k(end)' = 0
    
    g_end_mtx = [g_x(end); g_y(end)];
    k_end_mtx = [real(k(end)); imag(k(end))];
    g_proj = (k_end_mtx * k_end_mtx')/(k_end_mtx'*k_end_mtx)*g_end_mtx;
    g_proj_mag = sqrt(g_proj'*g_proj);  %target reduce this to 0
    
    blip_gr_sr_max = 9000; %G/cm/s  90 mT/m/ms
    blip_gr_sr_actual = -blip_gr_sr_max;
    gr_proj_str_step = blip_gr_sr_actual * gr_dwell;
    
    g_end_current_mtx = g_end_mtx;
    k_end_current_mtx = k_end_mtx;
    k_cpx = k_end_current_mtx(1) + 1i*k_end_current_mtx(2);
    finish = false;  step = 1;
    clear g_tailor_mtx k_tailor_mtx
    while(~finish)
        %=========visulize===================z
        %         figure(11);
        %
        %         quiver(0,0,k_end_current_mtx(1),k_end_current_mtx(2),'r')
        %         hold on
        %         quiver(k_end_current_mtx(1),k_end_current_mtx(2),g_end_current_mtx(1),g_end_current_mtx(2))
        %
        %         xlim([-2 2])
        %         ylim([-2 2])
        %         hold off
        %         drawnow(); pause(0.5);
        %====================================
        
        g_tailor_mtx(:, step) = g_end_current_mtx + gr_proj_str_step.*[cos(angle(k_cpx)); sin(angle(k_cpx))];
        k_tailor_mtx(:, step) = k_end_current_mtx + g_end_current_mtx.*gamma.* gr_dwell;
        
        g_end_current_mtx = g_tailor_mtx(:, step);
        k_end_current_mtx = k_tailor_mtx(:, step);
        k_cpx = k_end_current_mtx(1) + 1i*k_end_current_mtx(2);
        
        g_proj_current = (k_end_current_mtx * k_end_current_mtx')/(k_end_current_mtx'*k_end_current_mtx)*g_end_current_mtx;
        g_proj_current_mag = sqrt(g_proj_current'*g_proj_current);
        %         disp(['step: ',num2str(step),'; gr projected along k(end): ', num2str(g_proj_current_mag)]);
        
        if(g_end_current_mtx'*k_end_current_mtx/abs(k_cpx)<gr_proj_str_step/2)
            finish = true;
        else
            step = step +1;
        end
    end
    %set the last GR perpendicular to k
    g_tailor_mtx(:, step) = g_tailor_mtx(:, step) - g_proj_current;
    disp(['Smoothing tailor length: ', num2str(step)]);
    
    %     figure(12);
    %     subplot(131); plot(k_tailor_mtx(1,:),k_tailor_mtx(2,:)); title('trajecotry'); xlim([-1 1]); ylim([-1 1])
    %     subplot(132); plot(g_tailor_mtx(1,:)); hold on; plot(g_tailor_mtx(2,:),'r'); title('gradients');
    %     subplot(133); plot(abs(diff(g_tailor_mtx(1,:))./gr_dwell + 1i*diff(g_tailor_mtx(2,:))./gr_dwell)); title('SL')
    
    k_tailor = k_tailor_mtx(1,:) + 1i*k_tailor_mtx(2,:);
    k_optimized = [k k_tailor];
    
    phi_end = angle(k_optimized(end));
    k_optimized =  k_optimized .* exp(-i.*phi_end); %always ends at 0 degree
    
elseif(brutal)
    %%
    %initial
    phi_end = angle(k(end));
    k =  k .* exp(-i.*phi_end); %always ends at 0 degree
    g_y = 1/gamma * cat(2, 0, diff(imag(k))./gr_dwell);
    g_x = 1/gamma * cat(2, 0, diff(real(k))./gr_dwell);
    
    g_end_mtx = [g_x(end); g_y(end)];
    k_end_mtx = [real(k(end)); imag(k(end))];
    
    %max step
    blip_gr_sr_max = 9000; %G/cm/s  90 mT/m/ms
    blip_gr_sr_actual = -blip_gr_sr_max;
    gr_proj_str_step = blip_gr_sr_actual * gr_dwell;
    
    finish = false; step = 0;
    while(~finish)
        
        if(g_x(end)+gr_proj_str_step<0)
%             gr_proj_str_step = -g_x(end);
            finish = true;
        end
        
        g_end_mtx = [g_x(end)+gr_proj_str_step; g_y(end)];
        
        
        k_tailor_mtx = k_end_mtx + g_end_mtx.*gamma.* gr_dwell;
        k_tailor = k_tailor_mtx(1,:) + 1i*k_tailor_mtx(2,:);
        k = [k k_tailor];
        
        %update k
        phi_end = angle(k(end));
        k =  k .* exp(-i.*phi_end); %always ends at 0 degree
        k_end_mtx = [real(k(end)); imag(k(end))];
        figure(18); plot(k); drawnow;
        
        %update GR
        g_y = 1/gamma * cat(2, 0, diff(imag(k))./gr_dwell);
        g_x = 1/gamma * cat(2, 0, diff(real(k))./gr_dwell);
        [g_x(end); g_y(end)]
        step = step + 1;
        if(finish)
            disp(['steps = ', num2str(step)]);
        end
        
    end
    
    k_optimized = k;
    
    
end


end

