function  A = nuFTOperator(trajectory, imageDim, sensmaps, os)

% usage:
%    A = nuFTOperator(trajectory, imageDim, sensmaps)
%
% imageDim = [#rows #columns] = [Nr Nc], of the image
% trajectory = Nt x 2 = arbitrary trajectory in 2D k-space
% densityComp = variable density compensation vector (=1 if no comp.)
% sensmaps = Nr x Nc x #coils
% os = oversampling (default=1)

if nargin==0 % default constructor
    s.dim = [];
    s.numCoils = [];
    s.imageDim = [];
    s.adjoint = 0;
    s.trajectory_length = [];
    s.nufftNeighbors = [];
    s.sensmaps = {};
    s.nufftStruct = [];
else
    
    if size(trajectory,2) == 3 && length(imageDim)==3
        if imageDim(3) == 1
            s.dim = 2;
            trajectory = trajectory(:,1:2); %remove 3rd dimensions
            imageDim(3) = [];
        else
            %normal 3d case
            s.dim = 3;
        end
    elseif size(trajectory,2) == 2 && length(imageDim)==3
        %this makes only sense if imageDim(3) == 1
        if imageDim(3) == 1
            s.dim = 2;
            imageDim(3) = [];
        else
            warning('image dim and trajectory dim do not match.')
        end
    elseif size(trajectory,2) == 2 && length(imageDim)==2
        %normal 2d case
        s.dim = 2;
    end
     
    %check coil sensitivity dimension
    if nargin<=2 || isempty(sensmaps)
        s.numCoils = 1;
    else
        s.numCoils = size(sensmaps,s.dim+1);
    end
   
    %nufft oversampling
    if nargin<=3 || isempty(os)
        os = 1;
    end
    
    s.imageDim = imageDim;
    s.adjoint = 0;
    s.trajectory_length = size(trajectory,1);
    if size(trajectory,2) == 3
        s.nufftNeighbors = [5 5 5];
    elseif size(trajectory,2) == 2
        s.nufftNeighbors = [5 5];
    else
        s.nufftNeighbors = [5];
    end
    
    
    if nargin<=2 || isempty(sensmaps)
        s.sensmaps{1} = 1;
    else
        for k=1:s.numCoils
            if size(trajectory,2) == 3
                if s.numCoils > 1
                    s.sensmaps{k} = sensmaps(:,:,:,k);
                else
                    s.sensmaps{1}=sensmaps;
                end
            else
                s.sensmaps{k} = sensmaps(:,:,k);
            end
        end
    end
    
    if size(trajectory,2) == 3
       %trajectory = [trajectory(:,2), -trajectory(:,1) ,trajectory(:,3)]; %orginal gut und bewaehrt
       trajectory = [trajectory(:,1),  trajectory(:,2) ,trajectory(:,3)];
       n_shift = [imageDim(1)/2-1 imageDim(2)/2-1 round(imageDim(3)/2)-1]; %stimmt fuer ungerade Anzahl UND gerade Anzahl
    elseif size(trajectory,2) ==2
        %trajectory = [-trajectory(:,1), trajectory(:,2)];
        n_shift = [imageDim(1)/2-1 imageDim(2)/2-1];
    else
        trajectory = trajectory;
    end
    % the default fft shift is [N(1)/2 N(2)/2 N(3)/2]

    %the fft recon
    s.nufftStruct = nufft_init(trajectory, imageDim, s.nufftNeighbors, round(os*imageDim),n_shift,'kaiser');
%     s.nufftStruct = nufft_init(trajectory, imageDim, s.nufftNeighbors, round(os*imageDim));
end

A = class(s,'nuFTOperator');
