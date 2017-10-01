spatial_images=reshape(G,[sdu(1) sdu(2) size(G,2)]);

I1=[];
I2=[];
for i=1:size(G,2);
I1=    cat(2,I1,abs(spatial_images(:,:,i)));
I2=    cat(2,I2,angle(spatial_images(:,:,i)));
end

figure(9999);
subplot(211)
imshow(I1,[])
title('abs spatial images')

subplot(212)
imshow(I2,[-pi pi])
title('ph spatial images')
