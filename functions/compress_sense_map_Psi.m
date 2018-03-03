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
%
% (c) q.zhang @ AMC 2017

function [compressed_sense_map, compressed_sense_Psi] = compress_sense_map_Psi(VirtualCoilMartix, full_sense_map, full_sense_Psi)

assert(size(full_sense_map, 4) == size(full_sense_Psi, 1),'sense map does not match Psi');

full_sense_map_copy = full_sense_map;
cur_nr_chans = size(full_sense_map, 4);
for ac_in = 1:size( VirtualCoilMartix, 3 )
    row_ind = min( [size(VirtualCoilMartix, 1), find( isnan( VirtualCoilMartix(:,1,ac_in ) ),1 )-1]);
    col_ind = min( [size(VirtualCoilMartix, 2), find( isnan( VirtualCoilMartix(1,:,ac_in ) ),1 )-1]);
    A = VirtualCoilMartix(1:row_ind, 1:col_ind,ac_in);
    
    assert(size(A,2) == cur_nr_chans)
    
    full_sense_map = permute( full_sense_map, [4,1,2,3,5,6,7,8,9,10,11,12,13] );
    %     coil_ref = permute( coil_ref, [4,1,2,3,5,6,7,8,9,10,11,12,13] );
    siz = size(full_sense_map); siz(1) = size(A,1);
    compressed_sense_map = A*full_sense_map(:,:);
    %     coil_ref = A*coil_ref(:,:);
    compressed_sense_map = reshape(compressed_sense_map, siz);
    %     coil_ref = reshape(coil_ref, siz);
    compressed_sense_map = permute( compressed_sense_map, [2,3,4,1,5,6,7,8,9,10,11,12,13] );
    %     coil_ref = permute( coil_ref, [2,3,4,1,5,6,7,8,9,10,11,12,13] );
    compressed_sense_Psi = A*full_sense_Psi*A';
end


disp(['SENSE maps and Psi are compressed to ', num2str(size(A,1)) ,' virtual channels!']);
figure(102);
slice_id = max(1, floor(size(full_sense_map,3)/2));
subplot(121); montage(permute(squeeze(abs((full_sense_map_copy(:,:,slice_id,:)))),[1 2 4 3]),'displayrange',[]); title('SENSE before cc')
subplot(122); montage(permute(squeeze(abs((compressed_sense_map(:,:,slice_id,:)))),[1 2 4 3]),'displayrange',[]); title('SENSE after cc')
end