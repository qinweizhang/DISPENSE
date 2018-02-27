function image_filterd =  filt_perifiral_kspa(image_corrected_unfiltered, mask, trim_kspa_filter_mask_size)
%
%INPUT
%
%image_corrected_unfiltered: to-be-filtered image 
%mask:                       sampling mask coresponding to image_corrected_unfiltered
%
%(c) q.zhang 2017 @ AMC
%


[x_dim, y_dim,~]=size(image_corrected_unfiltered); %get first 2 dim
% assert(prod(size(mask)==size(image_corrected_unfiltered))==1,'mask should has the same size as image');  %original
assert(prod(size(mask)==[x_dim, y_dim])==1,'mask should has the same size as image'); 


mask_all = 0 * mask;
for x = 1:x_dim
    mask_idx = find(mask(x,:));
    if(~isempty(mask_idx))
      mask_all(x, mask_idx(1):mask_idx(end)) = 1;
    end
end

for y = 1:y_dim
    mask_idx = find(mask_all(:,y));
    if(~isempty(mask_idx))
      mask_all(mask_idx(1):mask_idx(end), y) = 1;
    end
end


%trim mask
mask_trim = trim_kspa_filter_mask_size; %possitive. Pixel numbers of the frame being discard.
if(mask_trim>0)
    mask_all = conv2(mask_all, ones(mask_trim,mask_trim)/(mask_trim.*mask_trim),'same');
    mask_all(find(mask_all<0.99)) = 0;
    mask_all(find(mask_all>0)) = 1;
end


%hard filter
kspace_window = mask_all;


%soft filter
beta = 100; %small=smoother
for x = 1:x_dim
    mask_idx = find(mask_all(x,:));
    if(~isempty(mask_idx))
        window_left = 0.5 + 1/pi*atan(beta*(([1:y_dim] - mask_idx(1))/round(abs(mask_idx(1)-y_dim/2))));
        window_right = 0.5 - 1/pi*atan(beta*(([1:y_dim] - mask_idx(end))/round(abs(mask_idx(1)-y_dim/2))));
        kspace_window(x, :) = window_left .* window_right;
    end
end

for y = 1:y_dim
    mask_idx = find(mask_all(:,y));
    if(~isempty(mask_idx))
        window_up = 0.5 + 1/pi*atan(beta*(([1:x_dim] - mask_idx(1))/round(abs(mask_idx(1)-x_dim/2))));
        window_down = 0.5 - 1/pi*atan(beta*(([1:x_dim] - mask_idx(end))/round(abs(mask_idx(1)-x_dim/2))));
        kspace_window(:, y) = kspace_window(:, y) .* window_up' .* window_down';
    end
end

figure(2001); 
subplot(121); imshow(mask,[]);title('sampling mask');
subplot(122); imshow(kspace_window,[]);title('kspace filter'); colormap jet
% 
% figure(2002);
% [X, Y] = meshgrid([1:y_dim],[1:x_dim]);
% mesh(X, Y, kspace_window)


kspace_filterd = bsxfun(@times, fft2d(image_corrected_unfiltered), kspace_window);
image_filterd = ifft2d(kspace_filterd);

