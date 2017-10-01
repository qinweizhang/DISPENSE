thickness = 54;
[ky,kz] = meshgrid(-54:54, -14:14);

% replace this with your deformation
% kx = cos(ky).*cos(kz) / 5;
kx = zeros(size(ky));


% c = kx;  %undersampling pattern
c = double(pp'>0);  %undersampling pattern


%===============plot


figure(500);

% top & bottom faces
surface(kx+thickness,ky,kz,c,'EdgeColor','none')
surface(kx-thickness,ky,kz,c,'EdgeColor','none')


% Now the 4 sides
surface([kx(1,:)+thickness; kx(1,:)-thickness],...
        [ky(1,:); ky(1,:)],...
        [kz(1,:); kz(1,:)], ...
        [c(1,:); c(1,:)])

surface([kx(end,:)+thickness; kx(end,:)-thickness],...
    [ky(end,:); ky(end,:)],...
        [kz(end,:); kz(end,:)], ...        
        [c(end,:); c(end,:)])

surface([kx(:,1)+thickness, kx(:,1)-thickness],...
    [ky(:,1), ky(:,1)],...
        [kz(:,1), kz(:,1)], ...        
        [c(:,1), c(:,1)])

surface([kx(:,end)+thickness, kx(:,end)-thickness],...
     [ky(:,end), ky(:,end)],...
        [kz(:,end), kz(:,end)], ...       
        [c(:,end), c(:,end)])
    
    
az = 70;
el = 20;
view(az, el);    
axis equal
colormap gray
set(gca,'xtick',[],'ytick',[],'ztick',[])