function res = LRT_cell_2_mat(kspace)

%This function converts incompatable kspace data to one matrix by zeros
%filling shorter frames. 
%
%Now only used for SoSNav - TSE recon

tse_k = cat(4, kspace{:,2});
tes_k_2d = reshape(tse_k, size(tse_k,1)*size(tse_k,2)*size(tse_k,3),size(tse_k,4));


nav_k = cat(3, kspace{:,1});
nav_k_pad = zeros(size(tse_k,1)*size(tse_k,2), size(nav_k,2), size(nav_k, 3));
nav_k_pad(1:size(nav_k,1),:,:) = nav_k;
nav_k_pad_2d = reshape(nav_k_pad, size(nav_k_pad,1)*size(nav_k_pad,2), size(nav_k_pad, 3));

res = cat(2, nav_k_pad_2d, tes_k_2d);

end