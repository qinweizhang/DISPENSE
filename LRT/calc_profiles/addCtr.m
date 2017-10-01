function mask= addCtr(mask, ctrsize)
assert(ndims(mask)==2); 
[ky,kz]=size(mask);
kyc=floor((1+ky)/2);
kzc=floor((1+kz)/2);

mask(kyc-ctrsize:kyc+ctrsize,kzc-ctrsize:kzc+ctrsize)=ones(2*ctrsize+1,2*ctrsize+1); 

end