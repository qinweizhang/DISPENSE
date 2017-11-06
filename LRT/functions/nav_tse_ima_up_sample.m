function ima_high_res = nav_tse_ima_up_sample(ima_low_res,high_dim, mask)

% [5]low resolution K space (Squared)
low_dim =  size(ima_low_res);
kspa_low_res = fft2d(ima_low_res);

% [4] high resolution K space (Squared)
kspace_high_res_square = zeros(high_dim(2), high_dim(2), low_dim(3), low_dim(4));
high_dim_square = size(kspace_high_res_square);
kx_range = floor(high_dim_square(1)-low_dim(1))/2+1 : floor(high_dim_square(1)-low_dim(1))/2+low_dim(1);
ky_range = floor(high_dim_square(2)-low_dim(2))/2+1 : floor(high_dim_square(2)-low_dim(2))/2+low_dim(2);
kspace_high_res_square(kx_range, ky_range, :, :)  = kspa_low_res;

% [3] high resolution Image (Squared)
im_high_res_cropped = ifft2d(kspace_high_res_square);

% [2] high resolution Image
ima_high_res  = zeros(high_dim(1), high_dim(2), low_dim(3), low_dim(4));
Ix_range =  floor(high_dim(1)-high_dim(2))/2+1 : floor(high_dim(1)-high_dim(2))/2+high_dim(2);
ima_high_res(Ix_range,:,:) = im_high_res_cropped ;
    
% [1] masked high resolution Image  
ima_high_res = bsxfun(@times, ima_high_res, mask); 
    
end