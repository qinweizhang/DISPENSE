%%  add simulated coil sensitivity 
function [Iout,S_normalized]=addcoilsensitvity_to_simulated(I,ncoils,varargin)

[nx ny ndim1,ndim2]=size(I);
assert(nx==ny)

if nargin == 3
    disp('using external sens map...')
    S = squeeze(varargin{1});
    assert(ncoils==size(S, 3));
    S_k = bart('fft -i 3', S);
    
    bart_comd_resize = sprintf('resize -c 0 %d 1 %d 2 %d', nx, ny, ncoils);
    S_k_rs = bart(bart_comd_resize, S_k);
    S = bart('fft 3', S_k_rs);
else
    S=bart(['phantom -S',num2str(ncoils),' -x',num2str(nx)]);
    S=squeeze(S); 
end


S_norm=sqrt(sum(abs(S).^2,3));
S_normalized=bsxfun(@rdivide,S,S_norm);

Iout= bsxfun(@times,permute(I,[1 2 5 3 4]),conj(S_normalized));
% Iout= bsxfun(@times,permute(I,[1 2 5 3 4]),(S_normalized));


% Iout=permute(Iout,[1 2 5 3 4]);
end