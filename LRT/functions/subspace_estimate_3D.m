function [nav_estimate_1,nav_estimate_2,eigenvals_1,eigenvals_2]= subspace_estimate_3D(kspace,params)
% used for estimatinbg the subspace from full 3D kspace (before running LRT
% recon) 


fprintf('Estimating subspace from 3D kspace...\n')
assert(ndims(kspace)==6,'dims') % kx ky kz nc param1 

mask=(abs((kspace(floor(size(kspace,1)/2),:,:,1,:,:))))>0;
% >>>>>>>>>>>>>>>>>>>>RECON FROM HERE<<<<<<<<<<<<<<<<<<<<<<<<<<<
[Kx1,Ky1,Kx2,Ky2]=findSharedKpoints(permute(mask,[2 3 5 6 1 4]),params);
fprintf('# of shared ky-kz points of dim 1: %i \n',numel(Kx1))
fprintf('# of shared ky-kz points of dim 2: %i \n',numel(Kx2))

%%%%%%% 2: estimate subspaces: generalized for non-square shared k-points

%% param dim 1
nav_parameter_dim1=[];
for iter=1:length(Kx1)
    nav_parameter_dim1=cat(1,nav_parameter_dim1,(kspace(:,Kx1(iter),Ky1(iter),:,:,params.subspacedim1)));
end
[nx ny nz nc npar]=size(nav_parameter_dim1); 
nav_parameter_dim1=reshape(nav_parameter_dim1,[nx*nc,npar]);

[nav_estimate_1,eigenvals_1]= subspace_estimator_multicoil(squeeze(nav_parameter_dim1),params.L3);
fprintf('subspace 1 finished \n')
%% param dim 2
nav_parameter_dim2=[];
for iter=1:length(Kx1)
    nav_parameter_dim2=cat(1,nav_parameter_dim2,(kspace(:,Kx2(iter),Ky2(iter),:,params.subspacedim2,:)));
end
[nx ny nz nc npar]=size(nav_parameter_dim2); 
nav_parameter_dim2=reshape(nav_parameter_dim2,[nx*nc,npar]);

[nav_estimate_2,eigenvals_2]= subspace_estimator_multicoil(squeeze(nav_parameter_dim2),params.L4);
fprintf('subspace 2 finished \n')

end
