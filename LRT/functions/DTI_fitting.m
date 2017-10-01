% INPUT
%
% DTI_data:    4D image data in size of [x, y, z, dir]
% g:           gradient table in size of [dir, 3]
% b:           b values in size of [dir, 1]
%
% OUTPUT
% MD map
% FA map
% eigen vectors
%
%(c) q.zhang . 2017 . AMC

function [MD, FA, Deigvec] = DTI_fitting(DTI_data, g, b)
%% validate data

[x, y, z, dir] = size(DTI_data);
assert(size(g,1)==dir && size(g,2)==3,'gradient table is not in the correct size!')
%%

I=ones(3,3);

%generate G matrix
gx=[g(2:end,1)];
gy=[g(2:end,2)];
gz=[g(2:end,3)];
G=[gx.*gx,2*gx.*gy,2*gx.*gz,gy.*gy,2*gy.*gz,gz.*gz];
%generate image matrix
%generate B
B=DTI_data(:,:,:,2:end);
B=bsxfun(@rdivide, B, DTI_data(:,:,1,1));
B=log(B)/b;
B(find(isnan(B)))=0;   %remove NaN
B(find(isinf(B)))=0;   %remove Inf


GG=inv(G'*G)*G'; %projection /least square oparator: for G*D = B;    original model: [u v i]*D*[u v i]' = B; D = [Dxx Dxy Dxz; Dyx Dyy Dyz; Dzx Dzy Dzz];

%% calculate D 
D=zeros(x,y,z,6);
Deigval = zeros(x,y,z,3);
Deigvec = zeros(x,y,z,3,3);

for iz=1:z
    for iy=1:y
        for ix=1:x
            d = GG*col(B(ix, iy, iz, :));
            D(ix, iy, iz,:)= d;
            if(~isnan(d)&~isinf(d))
                d_mtx = [d(1) d(2) d(3); d(2) d(4) d(5); d(3) d(5) d(6)];
                Deigval(ix, iy, iz,:)=eig(d_mtx);
                [Deigvec(ix, iy, iz,:,:), tmp2]=eig(d_mtx);
            end           
            
        end
    end
end


%% calculate FA & MD

FA=sqrt(3.*((Deigval(:,:,:, 1)-Deigval(:,:,:,2)).^2+(Deigval(:,:,:,2)-Deigval(:,:,:,3)).^2+(Deigval(:,:,:,3)-Deigval(:,:,:,1)).^2)) ./ ...
    sqrt(2.*(Deigval(:,:,:,1).^2+Deigval(:,:,:,2).^2+Deigval(:,:,:,3).^2));

MD=mean(abs(Deigval), 4);


%% display
try
figure(71);
imagesc(FA,[0 1]); title('FA'); axis off; axis equal;colorbar

figure(72); 
imagesc(MD, [0 max(MD(:))]); title('MD'); colormap jet; axis off; axis equal; colorbar

figure(73);
subplot(131); Color1=Deigvec(:,:,:,1); imagesc(abs(Color1)); title('eigen vector #1'); axis off; axis equal
subplot(132); Color2=Deigvec(:,:,:,2); imagesc(abs(Color2)); title('eigen vector #2'); axis off; axis equal
subplot(133); Color3=Deigvec(:,:,:,3); imagesc(abs(Color3)); title('eigen vector #3'); axis off; axis equal
drawnow;
catch
    disp('image display failed...');
end
end