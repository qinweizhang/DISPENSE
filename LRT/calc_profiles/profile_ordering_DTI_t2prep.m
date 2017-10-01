function profile_order=profile_ordering(mask,varargin)
%input: mask size ky,kz,ndim1,ndim2; radialflag; linearflag; visualizeflag
% linearflag: 0: row/ 1:column
if nargin==1
    radialflag=1; visualizeflag=0;
elseif nargin==2
    radialflag=varargin{1}; linearflag=0; visualizeflag=0;
elseif nargin==3
    radialflag=varargin{1}; linearflag=varargin{2};
elseif nargin==4
    radialflag=varargin{1}; linearflag=varargin{2}; visualizeflag=varargin{3};
elseif nargin==5
    radialflag=varargin{1}; linearflag=varargin{2}; visualizeflag=varargin{3}; nr_shot = varargin{4};
end

[ky,kz,nDim1,nDim2]=size(mask);
nr_points=sum(sum(mask(:,:,1,1)));
assert(xor(radialflag, linearflag),'Choose either radial or linear!');


vec= @(x) x(:);
% idea: for every shot: measure  points that are close in ky,kz after each
% other: for all shots sort points as such
% inter-shot: permute such that inter shot ordering is random
profile_order = [];
for dim2=1:nDim2;
    m=mask(:,:,:,dim2);
    for dim1=1:nDim1;
        
        [row,col]=find(mask(:,:,dim1,dim2));
        row=row-floor((1+ky)/2);
        col=col-floor((1+kz)/2);
        
        ETL = length(row)/nr_shot;
        
        %         figure(10); plot(row(index_r),col(index_r))
        
        if(linearflag)
            kp(1,:,dim1)=row;
            kp(2,:,dim1)=col;
        end
        
        
        if(radialflag)
            radial_cor = col + 1i.*row;
            for s = 1:nr_shot
                order{s} = [];
            end
            ETL = length(radial_cor)/nr_shot;
            %determine shot n0
            for k = 1:length(radial_cor)
                shot_num = ceil((angle(radial_cor(k))+pi)/(2*pi)*nr_shot);
                order{shot_num} = cat(1, order{shot_num}, radial_cor(k));
            end
            %rank order within shot
            pool = [];
            for s = 1:nr_shot
                order{s} = sort(order{s});
                %tailprofile if it's longer
                if(length(order{s})>ETL)
                    pool =cat(1, pool, order{s}(ETL+1:end));
                    order{s}(ETL+1:end) = [];
                end
            end
            for s = 1:nr_shot
                %patch the shorter ones
                l = length(order{s});
                if(l<ETL)
                    order{s} =cat(1, order{s}, pool(1:ETL-l));
                    pool(1:ETL-l) = [];
                end
            end
            order_final = [];
            for s = 1:nr_shot
                order_final = cat(1, order_final, order{s});
            end
            kp(1,:,dim1)=imag(order_final);
            kp(2,:,dim1)=real(order_final);
            
            
            
        end
        
        
    end
    %for DTI-T2prep: both dim1 & dim2 are in dynamics
    profile_order = cat(3, profile_order, kp);
    
    
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
    
    if(exist('ETL','var'))
        figure(4);
        for frame = 1:size(profile_order,3);
            for k=1:size(profile_order,2)
                im(profile_order(1,k,frame)+floor((1+ky)/2), profile_order(2,k,frame)+floor((1+kz)/2),1,frame) = mod(k-1,ETL)+1;
            end
        end
        im = reshape(im, size(im,1),size(im,2),nDim1, nDim2);
        immontage4D(im, []);colormap(jet(128)); colorbar; title('#rf'); title('ETL profile')
        figure(5);
        montage(permute(im(:,:,:,1),[1 2 4 3]),'displayrange',[]);colormap jet; title('ETL profile for the first volume')
        
    end
    
    
    
    
end








end

