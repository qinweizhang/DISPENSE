clear all; clc; close all;
%==========const.=================
max_sample = 18000; %Set in the PPE
gamma = 42.57746778; %(1/ms/mT)
gr_dwell = 6.4e-3; %in ms
gr_max_str = 40; % in mT/m
gr_max_slewRate = 50; %in mT/m/ms

%=========Scan Parameters============
Par.FOV_xy = 0.210;%FOV in M/P in m
Par.FOV_z = 0.040;  %FOV in S in m
Par.VOX_xy = 0.010; %Voxel size in M/P in m
Par.VOX_z = 0.010;  %Voxel size in S in m
%============END=====================

kz_max = 2 * pi / Par.VOX_z;
kxy_max = 2 * pi / Par.VOX_xy;

fun_Rz = @(kz, kz_max)(1+2 * abs(kz) ./ kz_max); %undersample profile in z (User define)
fun_Rr = @(r, kxy_max)(1.5 + abs(r) ./ kxy_max);  %undersample profile in z (User define)
%% Kz and Rz(undersample profile in z)
delta_kz = 2 * pi / Par.FOV_z;
kz = cat(2, fliplr(-delta_kz:-delta_kz:-kz_max),0:delta_kz:kz_max);
Rz = fun_Rz(kz,kz_max);

%model K = A*theta(t)*exp(-i*theta(t))
t = [0:gr_dwell:100]';
theta = 2 * t;
A = 10 * sqrt( 2 -1 * t./t(end));
K = A .* theta.*exp(-i .* theta);
Kx = real(K);
Ky = imag(K);
Gx = cat(1, 0, diff(Kx)./gr_dwell)./gamma;
Gy = cat(1, 0, diff(Ky)./gr_dwell)./gamma;
Gx_sr =  cat(1, 0, diff(Gx)./gr_dwell)./gamma;
Gy_sr = cat(1, 0, diff(Gy)./gr_dwell)./gamma;

figure(5); subplot(311); plot(t, Kx); hold on; plot(t, Ky);
subplot(312); plot(t, Gx); hold on; plot(t, Gy);
subplot(313); plot(t, Gx_sr); hold on; plot(t, Gy_sr);