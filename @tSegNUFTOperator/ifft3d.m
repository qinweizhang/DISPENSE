function res = ifft3d(x)

fctr = size(x,1)*size(x,2)*size(x,3);

res = zeros(size(x));


res = sqrt(fctr)*fftshift(fftshift(fftshift(ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(x,1),2),3),[],1),[],2),[],3),1),2),3);



end