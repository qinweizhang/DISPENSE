clear all; clc; close all;
%==========const.=================
%Units are: s, G, cm
gamma       = 4257.6;
max_sample  = 18000;    %Set in the PPE
gr_dwell    = 6.4e-6;   %in s
gr_max      = 4;        %in G/cm (40 mT/m)
sr_max      = 10000;    %in G/cm/s (100 mT/m/ms)


%=========Scan Parameters============
Par.FOV_xy      = 21; %FOV in M/P in cm
Par.FOV_z       = 12;  %FOV in S in cm
Par.VOX_xy      = 1;  %Voxel size in M/P in cm
Par.VOX_z       = 1;  %Voxel size in S in cm
Par.N           = 1;  %Interleaves
%FOV(r) = Sum ( Fcoeff(k)*(r/rmax)^(k-1))
Par.Fxy_coeff   = [1 -0.5] 	% FOVxy decreases linearly from 1 FOV to 0.5FOV. (R: 1 -> 2)
Par.Fz_coeff    = [1 -0.5] % FOVxy decreases linearly from 1 FOV to 0.4FOV. (R: 1 -> 2.5)
%============END=====================

kz_max  = 1 / Par.VOX_z / 2;
kxy_max = 1 / Par.VOX_xy /2;

%% Kz and (undersample profile in z)
kz_all = calc_kz(kz_max, Par.FOV_z, Par.Fz_coeff);

%% Calulate Kxy, Gxy and SRxy
k_all = [];
for kz_idx = 1:length(kz_all)
    
    %calc spiral shape
    if(mod(kz_idx, 2) ==1 ) %update kxy_max only for spiral-out (odd)trjectories.
        Fz_coeff_now = 0;
        for coeff = 1:length(Par.Fz_coeff)
            Fz_coeff_now = Fz_coeff_now+(Par.Fz_coeff(coeff))*abs((kz_all(kz_idx)/kz_max))^(coeff-1);
        end
        Fcoeff = Par.Fxy_coeff .* Fz_coeff_now .* Par.FOV_xy;
        kxy_max_now = sqrt(kxy_max.^2 - kz_all(kz_idx).^2);
    end
    [k,g,s,time] = vds(sr_max,gr_max,gr_dwell,Par.N,Fcoeff,kxy_max_now);
    
    %connect with previous spirals
    if(mod(kz_idx, 2) == 0) %even spirals; spiral in        
        %reverse order and flip in y axis(imaginary)
        k = fliplr(conj(k));
        %match the ends
        end_theta = angle(k_all(1,end) + i*k_all(2,end));
        phi = end_theta - angle(k(1)); 
        k = k .* exp(i.*phi);    
        %connect the ends
        d_kz = kz_all(kz_idx) - kz_all(kz_idx - 1);
        blip_sr = 2000; %G/cm/s  20 mT/m/ms
        slop = sqrt(d_kz / gamma /blip_sr); %seconds
        gz = [0:gr_dwell:slop].*blip_sr;
        gz = [gz fliplr(gz(1:end-1))];
        kz_connect = k_all(3,end);
        for idx = 1:length(gz)
            kz_connect = [kz_connect kz_connect + gamma*gr_dwell*gz(idx)];
        end
        
        
    end
    k_3D = cat(1, real(k), imag(k), ones(1, length(k)).*kz_all(kz_idx) );
    k_all = [k_all k_3D];
    
    
end
figure(2); plot3(k_all(1,:),k_all(2,:),k_all(3,:));

% disp('Plotting Gradient');
% g_all = [real(g_all(:)) imag(g_all(:))];
% plotgradinfo(g_all,gr_dwell);


