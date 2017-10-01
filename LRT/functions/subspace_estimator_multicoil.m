function [nav_estimate,eigenvals]= subspace_estimator_multicoil(kspace,L)
% Function that calculates the subspace based on k-space
% estimator 
% to do: make compatible with 3D kspace 
%input:
%kspace: 3D kspace center size(kx,ky,nc,M)/(kx*ky,nc,M)
%L:  rank of subspace estimator (L) (L<=M)

if ndims(kspace)==4 % kx ky nc param
    assert(L<=size(kspace,4))
    S=reshape(kspace,[size(kspace,1)*size(kspace,2)*size(kspace,3),size(kspace,4)]);

elseif ndims(kspace)==3 %kx ky param
    kspace=permute(kspace,[1 2 4 3]);  
    assert(L<=size(kspace,4))
    S=reshape(kspace,[size(kspace,1)*size(kspace,2)*size(kspace,3),size(kspace,4)]);

elseif ndims(kspace)==2 %assume kx*ky*nx, param (already Casorati matrix)
    S=kspace; 
end  
    
% calculate singular value decomposition
[left_1,eigen_1,right_1]=svd(S.','econ');
eigenvals=diag(eigen_1);  %output for evaluation

%navigator estimate are first L left singular vectors
nav_estimate=left_1(:,1:L);



end
