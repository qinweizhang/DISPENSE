function  kspace_low_res = tse_nav_kspa_down_sample(kspace_high_res, low_dim)

% [1] high resolution k space
high_dim = size(kspace_high_res);

% [2] high resolution Image
im_high_res = ifft2d(kspace_high_res);

% [3] high resolution Image (Squared)
Ix_range =  floor(high_dim(1)-high_dim(2))/2+1 : floor(high_dim(1)-high_dim(2))/2+high_dim(2);
im_high_res_cropped = im_high_res(Ix_range,:,:);


% [4] high resolution K space (Squared)
kspace_high_res_square = fft2d(im_high_res_cropped);
high_dim_square =  size(kspace_high_res_square);

% [5]low resolution K space (Squared)
kx_range = floor(high_dim_square(1)-low_dim(1))/2+1 : floor(high_dim_square(1)-low_dim(1))/2+low_dim(1);
ky_range = floor(high_dim_square(2)-low_dim(2))/2+1 : floor(high_dim_square(2)-low_dim(2))/2+low_dim(2);
kspace_low_res = kspace_high_res_square(kx_range, ky_range, :);
    
end
