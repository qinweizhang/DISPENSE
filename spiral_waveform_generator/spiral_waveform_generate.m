clear all; clc; close all;
%==========const.=================
%Units are: s, G, cm
gamma       = 4257.6;
max_sample  = 18000;    %Set in the PPE
gr_dwell    = 6.4e-6;   %in s
gr_max      = 4;        %in G/cm (40 mT/m)
sr_max      = 10000;    %in G/cm/s (100 mT/m/ms)


%=========Scan Parameters============
Par.FOV_xy      = 22; %FOV in M/P in cm
Par.FOV_z       = 12;  %FOV in S in cm
Par.VOX_xy      = 0.5;  %Voxel size in M/P in cm
Par.VOX_z       = 0.5;  %Voxel size in S in cm
Par.N           = 1;  %Interleaves
%FOV(r) = Sum ( Fcoeff(k)*(r/rmax)^(k-1))
Par.Fxy_coeff   = [0.4 -0.2] 	% FOVxy decreases linearly from 1 FOV to 0.5FOV. (R: 1 -> 2)
Par.Fz_coeff    = [0.3] % FOVxy decreases linearly from 1 FOV to 0.4FOV. (R: 1 -> 2.5)
%============END=====================

kz_max  = 1 / Par.VOX_z / 2;
kxy_max = 1 / Par.VOX_xy /2;

%% Kz and (undersample profile in z)
kz_all = calc_kz(kz_max, Par.FOV_z, Par.Fz_coeff);
kz_all = -1*kz_all;
%% Calulate Kxy, Gxy and SRxy
k_all = [];
for kz_idx = 1:length(kz_all)
    %calc spiral shape
    if(mod(kz_idx, 2) ==1 ) %update kxy_max only for spiral-out (odd)trjectories.
        Fz_coeff_now = 0;
        for coeff = 1:length(Par.Fz_coeff)
            Fz_coeff_now = Fz_coeff_now+(Par.Fz_coeff(coeff))*abs((kz_all(kz_idx)/kz_max))^(coeff-1);
        end
        %         Fcoeff = Par.Fxy_coeff .* Fz_coeff_now .* Par.FOV_xy;  %original
        Fcoeff = Par.Fxy_coeff .* 1 .* Par.FOV_xy;
%         kxy_max_now = sqrt(kxy_max.^2 - kz_all(kz_idx).^2);
        kxy_max_now = kxy_max * sqrt(1 - (1-max(1-abs(kz_all(kz_idx)/kz_max),0.1)).^2) * kz_max;

    end
    [k,g,s,time] = vds(sr_max,gr_max,gr_dwell,Par.N,Fcoeff,kxy_max_now);

    assert(length(k)>0);
    k = smooth_traj_end(k);
    phi_end = angle(k(end)) + kz_idx.*pi/4;
    k =  k .* exp(-i.*phi_end); %always ends at 0 degree
    
    %connect with previous spirals
    if(mod(kz_idx, 2) == 0) %even spirals; spiral in
        %reverse order and flip in y axis(imaginary)
        k = fliplr(conj(k));
        %match the ends
        end_theta = angle(k_all(1,end) + i*k_all(2,end));
        phi = end_theta - angle(k(1));
        k = k .* exp(i.*phi);
        %connect the ends
        kz_connect = connect_kz(kz_all(kz_idx), kz_all(kz_idx - 1));
        assert((length(k)-length(kz_connect))>=0,'spiral-out at this kz location is not long enough for kz connection.');
        %                 k_3D = cat(1, real(k), imag(k), [ones(1, (length(k)-length(kz_connect)-round((length(k)-length(kz_connect))/2))).*kz_all(kz_idx-1), kz_connect, ones(1, round((length(k)-length(kz_connect))/2)).*kz_all(kz_idx)] );
%         k_3D = cat(1, real(k), imag(k), [ones(1, length(k)-length(kz_connect)).*kz_all(kz_idx-1) kz_connect] );
        k_3D = cat(1, real(k), imag(k), [kz_connect ones(1, length(k)-length(kz_connect)).*kz_all(kz_idx)] );
    else
        %connect the ends
        k_3D_connect = [];
        if(kz_idx > 1)
            kz_connect = connect_kz(kz_all(kz_idx), kz_all(kz_idx - 1));
            k_3D_connect = cat(1, zeros(2, length(kz_connect)), kz_connect);
        end
        k_3D = cat(1, real(k), imag(k), ones(1, length(k)).*kz_all(kz_idx) );
        k_3D = [k_3D_connect k_3D];
    end
    
    k_all = [k_all k_3D(:,2:end)];
    
end

%initial ramp down
kz_rampdown = connect_kz(kz_all(1), 0);
k_3D_rampdown = cat(1, zeros(2, length(kz_rampdown)), kz_rampdown);
k_all = [k_3D_rampdown k_all];

figure(2); plot3(k_all(1,:),k_all(2,:),k_all(3,:));
for dir = 1:3
    g_all(:,dir) = 1/gamma * cat(2, 0, diff(k_all(dir,:))./gr_dwell);
end
for dir = 1:3
    sr_all(:,dir) = cat(2, 0, diff(g_all(:,dir))'./gr_dwell);
end

disp('Plotting Gradient');
figure(10);
plotgradinfo(g_all(:,[1 2]),gr_dwell);

clc
disp(['samples: ' num2str(length(g_all))]);
disp(['total dur: ' num2str((length(g_all)-1)*gr_dwell*1000), ' ms']);
%% export


G_cm2mT_m = 10; %unit convert from G/cm to mT/m
G_cm_s2mT_m_ms = 1/100;%unit convert from G/cm/s to mT/m/ms

g_all_export = g_all * G_cm2mT_m;
if(ispc)    
    dlmwrite('L:\basic\divi\Ima\parrec\Kerry\External_spiral_waveform\external_spiral_GRwaveforms.dat',g_all_export,'delimiter','\t');
else
    dlmwrite('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/External_spiral_waveform/external_spiral_GRwaveforms.dat',g_all_export,'delimiter','\t');
end
