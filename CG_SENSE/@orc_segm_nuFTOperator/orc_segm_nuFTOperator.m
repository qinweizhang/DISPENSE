function  A = orc_segm_nuFTOperator(trajectory, imageDim, sensmaps, phasemap, dwelltime, Ns)

% function  A = orc_segm_nuFTOperator(trajectory, imageDim, sensmaps, phasemap, dwelltime, Ns)
%
% imageDim = [#rows #columns] = [Nr Nc], of the image
% trajectory = Nt x 2 = arbitrary trajectory in 2D k-space
% sensmaps = Nr x Nc x #coils
% phasemap = offresonance map in rads/s
% dwelltime = time between k-space samples
% Ns = number of segments


if isempty(sensmaps)
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

s.imageDim = imageDim;
s.adjoint = 0;
tl = size(trajectory,1);
s.trajectory_length = tl;
if size(trajectory,2) == 3
    s.nufftNeighbors = [5 5 5];
else
    s.nufftNeighbors = [5 5];
end
s.phasemap = phasemap;
s.tADC = tl*dwelltime;
s.oversampling = 2*s.imageDim;

s.nSegments = Ns;
sr = floor((tl-Ns-2)/Ns);
sl = 2*sr + 1;
sr_1 = sr;
sr_end = rem(tl-Ns-2,Ns);

si = cell(1,Ns+2);
si{1} = [1:sr_1+1];
for k=1:Ns-1
    si{k+1} = [k*(sr+1)+1-sr:k*(sr+1)+1+sr];
end
si{end-1} = [Ns*(sr+1)+1-sr:Ns*(sr+1)+1+sr_end];
si{end} = [tl-sr_end:tl];
s.segment_index = si;

h = hann(sl+2);
h = h(2:end-1);
h2 = hann(2*sr_end+1+2);
h2 = h2(2:end-1);
ipf = cell(1,Ns+2);
ipf{1} = h(sr+1:end);
for k=1:Ns-1
    ipf{k+1} = h;
end
ipf{end-1} = [h(1:sr);h2(sr_end+1:end)];
ipf{end} = h2(1:sr_end+1);
s.segment_filter = ipf;

T = s.tADC * [0:tl-1]/tl;
t = zeros(1,Ns+2);
t(1) = 0;
for k=1:Ns
    t(k+1) = T(k*(sr+1)+1);
end
t(end) = T(end);
s.t = t;

if isempty(sensmaps)
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
       n_shift = [imageDim(1)/2-1 imageDim(2)/2-1 round(imageDim(3)/2)-1]
    elseif size(trajectory,2) ==2
        trajectory = [trajectory(:,1), trajectory(:,2)];
        n_shift = [imageDim(1)/2-1 imageDim(2)/2-1];
    else
        trajectory = trajectory;
    end
    % the default fft shift is [N(1)/2 N(2)/2 N(3)/2]
    
    %oversampling
    os =2;
    
for k=1:Ns+2
    %nstr = nufft_init(trajectory(si{k},:), s.imageDim, s.nufftNeighbors, s.oversampling, s.imageDim/2, 'kaiser');
    nstr = nufft_init(trajectory(si{k},:), s.imageDim, s.nufftNeighbors, s.oversampling,n_shift,'kaiser');
    s.interpolation_matrix{k} = nstr.p.arg.G;
end
s.scaling_factor = nstr.sn; % is the same each time

A = class(s,'orc_segm_nuFTOperator');
