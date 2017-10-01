function profile_order=profile_ordering(mask,varargin)
%input: mask size ky,kz,ndim1,ndim2; radialflag; linearflag; visualizeflag
% linearflag: 0: row/ 1:column
if nargin==1
    radialflag=0; visualizeflag=0;
elseif nargin==2
    radialflag=varargin{1}; linearflag=0; visualizeflag=0;
elseif nargin==3
    radialflag=varargin{1}; linearflag=varargin{2};
else
    radialflag=varargin{1}; linearflag=varargin{2}; visualizeflag=varargin{3};
end

[ky,kz,nDim1,nDim2]=size(mask);
nr_points=sum(sum(mask(:,:,1,1)));



vec= @(x) x(:); 
% idea: for every shot: measure  points that are close in ky,kz after each
% other: for all shots sort points as such 
% inter-shot: permute such that inter shot ordering is random 

for dim2=1:nDim2;
    m=mask(:,:,:,dim2);
    for dim1=1:nDim1;
        [row,col]=find(mask(:,:,dim1,dim2));
        row=row-floor((1+ky)/2);
        col=col-floor((1+kz)/2);

        % sort into a certain number of radial spokes: low to high 
        [theta, rho]=cart2pol(row,col);
        
        if radialflag
        [~,s_index_angle]=sort(theta);
        else
            if linearflag==0
                [~,s_index_angle]=sort(row);
            else
                [~,s_index_angle]=sort(col);
            end
            
        end

        % take blocks of N angle pieces (plus a residual) 
        rho_sorted=col(s_index_angle); 
        
        N=15;
        number_of_pieces=floor(nr_points/N);
        residual_N=mod(nr_points,N);
        
        index_r=[];
        for ii=1:number_of_pieces;
        [~,index_piece]=sort(rho_sorted(1+(ii-1)*N:ii*N));
        index_piece_transform=s_index_angle((ii-1)*N+index_piece);
        index_r=[index_r;index_piece_transform];
        end
        [~,index_piece]=sort(rho_sorted(1+number_of_pieces*N:end));
        index_piece_transform=s_index_angle((number_of_pieces*N)+index_piece);
        index_r=[index_r;index_piece_transform];

%         figure(10); plot(row(index_r),col(index_r))

        kp(1,:,dim1)=row(index_r);
        kp(2,:,dim1)=col(index_r);
        
        % sort every block of alpha cols form low to high in rows
        
        
    end
    
    profile_order(1,:,dim2)=vec(permute(kp(1,:,:),[1 3 2]));
    profile_order(2,:,dim2)=vec(permute(kp(2,:,:),[1 3 2]));
    
end



%% visualize ordering 
if visualizeflag
figure(2);clf
 %plot 100 TSE trains
for i=1:nDim1:nDim1*100; 
    hold on 
    plot(profile_order(1,1:i),profile_order(2,1:i),'k.')
    plot(profile_order(1,i:i+60),profile_order(2,i:i+60),'r.')
    hold off
    xlim([-64 64])   
    ylim([-64 64])
    pause(0.1)
    drawnow;
end 

%% 
clear delta_x delta_y
p=profile_order;
delta_x=(p(1,1:end-1,1)-p(1,2:end,1));
delta_y= (p(2,1:end-1,1)-p(2,2:end,1));

figure(3); clf; 
subplot(211);
plot(sqrt((delta_x).^2+(delta_y).^2),'.');
subplot(212);
plot((delta_x),'r*'); hold on
plot((delta_y),'b*'); hold off
title('kspace-jumps')





end








end

