function res = fft3d(x)
fctr = size(x,1)*size(x,2)*size(x,3);
res = zeros(size(x));


res = (1/sqrt(fctr))*fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(x,1),2),3),[],1),[],2),[],3),1),2),3);

% res = reshape(res, size_x);

end