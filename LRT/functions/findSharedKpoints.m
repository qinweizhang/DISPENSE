function [Kx1,Ky1,Kx2,Ky2]=findSharedKpoints(mask,params)
%input: mask (kx,ky,dims1,dim2)
% goal: find points which are sampled in all dim1,1 and 1,dim2

Cdim1=squeeze(sum(mask(:,:,:,params.subspacedim1),3))==size(mask,3);
assert(sum(Cdim1(:))>0,'no shared k-space points found for parameter dim 1!')
Cdim2=squeeze(sum(mask(:,:,params.subspacedim2,:),4))==size(mask,4);
assert(sum(Cdim2(:))>0,'no shared k-space points found for parameter dim 2!')

[Kx1,Ky1]=find(Cdim1);
[Kx2,Ky2]=find(Cdim2);
end