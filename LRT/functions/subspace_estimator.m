function nav_estimate= subspace_estimator(kspace,L)
% Function that calculates the subspace based on k-space
% estimator 
% to do: make compatible with 3D kspace 
disp('estimating subspace...')
%input:
%kspace: 3D kspace center size(kx,ky,M)
%L:  rank of subspace estimator (L) (L<=M)
assert(L<=size(kspace,3))

%reshape to Casorati matrix
S=reshape(kspace,[size(kspace,1)*size(kspace,2),size(kspace,3)]);

% calculate singular value decomposition
[left_1,eigen_1,right_1]=svd(S);

%navigator estimate are first L left singular vectors
nav_estimate=right_1(:,1:L);

end
