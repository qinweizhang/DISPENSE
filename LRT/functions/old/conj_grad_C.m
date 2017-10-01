function Ck=conj_grad_C(Gk,C,Ak,Bk,Y,Z,beta,du,Phi)
%Lambda: Omega with var C

BZ=(beta/2).*(Bk+Z./beta);

zf=fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(du,1),2),[],1),[],2),1),2); %zero filled recon
zf=reshape(zf,[size(zf,1)*size(zf,2),size(zf,3)*size(zf,4)]);
%zf = G C Phi
C_hat=Gk'*zf*(Phi);

numerator=C_hat+BZ;
denom=(1+beta/2).*ones(size(C_hat)); % TO DO?/
denom=denom.^(-1);
disp('únsure how to deal with Lambda....' )

Ck=denom.*numerator;
end