function  A = nuFTOperator(trajectory, imageDim, sensmaps,saveMemory)

% usage:
%    A = nuFTOperator(trajectory, imageDim, sensmaps)
%
% imageDim = [#rows #columns] = [Nr Nc], of the image
% trajectory = Nt x 2 = arbitrary trajectory in 2D k-space
% densityComp = variable density compensation vector (=1 if no comp.)
% sensmaps = Nr x Nc x #coils
% os = oversampling (default=2)

%still Thimo Hugger's (formerly known as Grotz) nufft-operator for CG-SENSE
%adapted for large 3D retrospective motion correction problems


if nargin==0 % default constructor
    s.numCoils = [];
    s.imageDim = [];
    s.adjoint = 0;
    s.trajectory_length = [];
    s.nufftNeighbors = [];
    s.sensmaps = {};
    s.nufftStruct = [];
else
    if nargin<=2 || isempty(sensmaps)
        s.numCoils = 1;
    else
        if size(trajectory,2) == 3 && length(size(sensmaps))== 4
            s.numCoils = size(sensmaps, length(size(sensmaps)));
        end
        if size(trajectory,2) == 3 && length(size(sensmaps))== 3
            s.numCoils = 1;
        end
        if size(trajectory,2) == 2 && length(size(sensmaps))== 3
            s.numCoils = size(sensmaps, length(size(sensmaps)));
        end
        if size(trajectory,2) == 2 && length(size(sensmaps))== 2
            s.numCoils = 1;
        end
    end
    if nargin<=3
        saveMemory = false;
    end
    
    os = 2; %we always want to have two-fold oversmapling
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
       trajectory = [trajectory(:,1),  trajectory(:,2) ,trajectory(:,3)];
       n_shift = [imageDim(1)/2-1 imageDim(2)/2-1 round(imageDim(3)/2)-1]; %stimmt fuer ungerade Anzahl UND gerade Anzahl
    elseif size(trajectory,2) ==2
        %trajectory = [-trajectory(:,1), trajectory(:,2)];
        n_shift = [imageDim(1)/2-1 imageDim(2)/2-1];
    else
        trajectory = trajectory;
    end
    % the default fft shift is [N(1)/2 N(2)/2 N(3)/2]

    %traj=reshape(trajectory,[nd*nl*ns*nt*3])*2*pi;
   
        %the fft recon
        if saveMemory
            %use table lookup for interpolation kernels
            s.nufftStruct = nufft_init(trajectory, imageDim, s.nufftNeighbors, round(os*imageDim),n_shift,'table',2^20,'minmax:kb');
            display('Using memory efficient table lookup for nufft.')
        else
            %use precomputed kernel
            s.nufftStruct = nufft_init(trajectory, imageDim, s.nufftNeighbors, round(os*imageDim),n_shift);
        end
        
end

A = class(s,'nuFTOperator');
