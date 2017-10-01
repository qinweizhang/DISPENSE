% recon simulation
% draft reconstruction code 
clear all; close all; clc;

% 1: make data (most settings in other .m file for now)
uf=0.12; % undersampling factor (excluding center)
noiselevel=0.3;
ncoils=5;
complexsim=1;
center_for_all_frames=0;
load('sense_map.mat'); %optional load external sense map

simulation_typ = 'qMATCH';
run create_undersampled_measurement_S.m    


%% SENSE map estimation
% a = sum(sum(du(:,:,:,:,:),4),5)./sum(sum(du(:,:,:,:,:)~=0,4),5);
% a =du(:,:,:,1,1);
a =d(:,:,:,1,1);
sens_est=bart('ecalib -m1 -c0',permute(a,[4 1 2 3]));
% sens_est=sens_est+1e-7; % no zero vals in sense maps...

% sens_est = bart('caldir 10',permute(a,[4 1 2 3]));
% sens_est=sens;
figure(112)
subplot(321)
immontage4D(abs(sens_est),[]); title('est. sense abs')
subplot(323)
immontage4D(real(sens_est),[]); title('est. sense real')
subplot(325)
immontage4D(angle(sens_est),[-pi pi]); title('est. sense phase')
subplot(322)
immontage4D(permute(abs(sens),[4 1 2 3]),[]); title('true sense abs')
subplot(324)
immontage4D(permute(real(sens), [4 1 2 3]),[]); title('true sense real')
subplot(326)
immontage4D(permute(angle(sens), [4 1 2 3]),[-pi pi]); title('true sense phase')

figure(113)
subplot(311)
immontage4D(permute(abs(squeeze(sens_est)./sens),[4 1 2 3]),[0.9 1.1]); title('true sense abs')
subplot(312)
immontage4D(permute(real(squeeze(sens_est)./sens), [4 1 2 3]),[0.9 1.1]); title('true sense real')
subplot(313)
immontage4D(permute(angle(squeeze(sens_est)./sens), [4 1 2 3]),[-pi pi]); title('true sense phase'); colormap jet

%% 
close all;
params=params_init();
params.Lg=2;
params.L3=5;
params.L4=4;
params.sparsity_transform='TV';
params.Imref=I;
params.x=95;
params.y=60;
params.increase_penalty_parameters=false;
params.inspectLg=true;
params.subspacedim1=1;
params.subspacedim2=1; 
params.G.precon=true;
params.G.maxiter = 10;
params.scaleksp=1;

% params.mu=1e3;
% params.lambda=5e-4;
% P_recon=LRT_recon(du,squeeze(sens_est),params);
P_recon=LRT_recon(du,squeeze(sens),params);