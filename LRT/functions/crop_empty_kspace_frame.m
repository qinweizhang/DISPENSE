function kspace_corpped = crop_empty_kspace_frame(kspace)

 filter_mask = sum(sum(squeeze(abs(kspace(:,:,1,:,:)))>0,3),4)>0;


x_status = sum(filter_mask, 1)>0;
x_valid = find(x_status);

y_status = sum(filter_mask, 2)>0;
y_valid = find(y_status);

kspace_corpped = kspace(y_valid(1):y_valid(end),x_valid(1):x_valid(end),:,:,:);

end
