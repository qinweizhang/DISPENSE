% Orthogonal SENSE maps: recombine coils to make the Psi map indentity mtx (SNR optimized)
% 
% 
% (c) q.zhang

function [kspa_orthocoil, sense_map_orthocoil] = Ortho_data_sense(sense_Psi, kspa, sense_map)


for t=1:length(sense_Psi)  %make sure diag(sense_Psi) are all real; they should be!
    sense_Psi(t,t) = real(sense_Psi(t,t));
end
L = chol(sense_Psi,'lower'); %Cholesky decomposition; L: lower triangle
L_inv = inv(L);
for c = 1:size(sense_Psi,1)
    %recombine sense map
    sense_map_orthocoil(:,:,c) = sum(bsxfun(@times, sense_map, permute(conj(L_inv(c,:)),[1 3 2])), 3);
    %recombine kspa map
    kspa_orthocoil(:,:,c,:) = sum(bsxfun(@times, kspa, permute(conj(L_inv(c,:)),[1 3 2])), 3);
end
figure(401);
subplot(121); montage(permute(abs(sense_map),[1 2 4 3]),'displayrange',[]); title('originial SENSE map')
subplot(122); montage(permute(abs(sense_map_orthocoil),[1 2 4 3]),'displayrange',[]); title('Orthogonal SENSE map')

%renormalize sense
sense_map_orthocoil = squeeze(normalize_sense_map(permute(sense_map_orthocoil,[1 2 4 3])));

