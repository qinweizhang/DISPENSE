% INPUT
% 
% VirtualCoilMartix: [N_virtual_c, N_full_c]
% full_sense_map:    [Nx, Ny, Nz, N_full_c]
% full_sense_Psi:    [N_full_c, N_full_c]; noise covarance mtx
% 
% OUT
% 
% compressed_sense_map:  [Nx, Ny, Nz, N_virtual_c]
% compressed_nav_sense_Psi: [N_virtual_c, N_virtual_c]

function [compressed_sense_map, compressed_nav_sense_Psi] = compress_sense_map_Psi(VirtualCoilMartix, full_sense_map, full_sense_Psi)

assert(size(full_sense_map, 4) == size(full_sense_Psi, 1),'sense map does not match Psi');
assert(size(full_sense_map, 4) == size(VirtualCoilMartix, 2),'sense map does not match Convert Matrix');

[N_virtual_c, ~] = size(VirtualCoilMartix);


compressed_sense_map = zeros(size(full_sense_map, 1),size(full_sense_map, 2),size(full_sense_map, 3),N_virtual_c);
for c = 1:N_virtual_c
    compressed_sense_map(:,:,:,c) = sum(bsxfun(@times, full_sense_map, permute(VirtualCoilMartix(c,:),[4 3 1 2]) ) , 4);
end

compressed_nav_sense_Psi =  VirtualCoilMartix * full_sense_Psi * VirtualCoilMartix'; %?  or VirtualCoilMartix.'



end
