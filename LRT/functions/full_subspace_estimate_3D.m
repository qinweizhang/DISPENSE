function [nav_estimate_1,eigenvals_1,nav_estimate_2,eigenvals_2] = full_subspace_estimate_3D(kspace,params)

% needed: params.L3 /L4
% kspace (kx,ky,kx,coils,params1,params2)

% find mask & shared k-space points 
mask=squeeze(kspace(1,:,:,1,:,:))~=0;
[Kx1,Ky1,Kx2,Ky2]=findSharedKpoints(mask,params);
fprintf('# of shared ky-kz points of dim 1: %i \n',numel(Kx1))
fprintf('# of shared ky-kz points of dim 2: %i \n',numel(Kx2))

figure(1); clf; immontage4D(mask,[0 1]);
xlabel('Parameter 1'); ylabel('Parameter 2');


% 2: estimate subspaces: generalized for non-square shared k-points
nav_parameter_dim1=[];
for iter=1:length(Kx1)
    nav_parameter_dim1=cat(1,nav_parameter_dim1,(kspace(Ky1(iter),Ky1(iter),:,:,params.subspacedim1)));
end
[nav_estimate_1,eigenvals_1]= subspace_estimator_multicoil(squeeze(nav_parameter_dim1),params.L3);

nav_parameter_dim2=[];
for iter=1:length(Kx2)
    nav_parameter_dim2=cat(1,nav_parameter_dim2,(kspace(Kx2(iter),Ky2(iter),:,params.subspacedim2,:)));
end
[nav_estimate_2,eigenvals_2]= subspace_estimator_multicoil(squeeze(nav_parameter_dim2),params.L4);


