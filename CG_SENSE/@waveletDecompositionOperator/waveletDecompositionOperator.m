function  W = waveletDecompositionOperator(dim, N, wname)

% function W = waveletDecompositionOperator(dim, N, wname)
%
% dim = dimension of image (e.g. [64 64 32])
% N = order of decomposition
% wname = name of wavelet

if nargin==0
    s.space = 0;
else
    if length(dim)==1
        s.space = 1;
    elseif length(dim)==2
        if dim(1)==1 || dim(2)==1
            s.space = 1;
        else
            s.space = 2;
        end
    elseif length(dim)==3
        s.space = 3;
    end
end
        
if s.space==0
    s.params.sizes = [];
    s.params.sizeINI = [];
    s.params.level = [];
    s.params.mode = [];
    s.params.filters = [];
    s.dim = [];
    s.space = [];
    s.sizes_all = [];    
elseif s.space==1
    [~,sizes] = wavedecX(zeros(dim,1),N,wname);
    s.params.sizes = sizes;
    s.params.sizeINI = [];
    s.params.level = [];
    s.params.mode = [];
    s.params.filters = [];
    s.dim = [dim 1];
    s.space = 1;
    s.sizes_all = [];
elseif s.space==2
    [~,sizes] = wavedec2X(zeros(dim),N,wname);
    s.params.sizes = sizes;
    s.params.sizeINI = [];
    s.params.level = [];
    s.params.mode = [];
    s.params.filters = [];
    s.dim = dim;
    s.space = 2;
    s.sizes_all = [];
elseif s.space==3
    s.params = wavedec3X(zeros(dim),N,wname);
    s.dim = dim;
    s.space = 3;
    sizes_all = zeros(length(s.params.dec),3);
    for k=1:length(s.params.dec)
        sizes_all(k,:) = size(s.params.dec{k});
    end
    s.sizes_all = sizes_all;
else
    error('dimension is not specified properly.');
end

if nargin==0
    s.N = [];
    s.wname = '';
    s.adjoint = 0;    
else
    s.N = N;
    s.wname = wname;
    s.adjoint = 0;
end

W = class(s,'waveletDecompositionOperator');