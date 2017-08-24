function Q = mtimes_par(A,B)

tra_length = A.trajectory_length;

if strcmp(class(A),'nuFTOperator')
    if A.adjoint
        Q = zeros(size(A));
        Btemp = reshape(B,[tra_length A.numCoils]);
        Anufft = A.nufftStruct;
        sn = single(Anufft.sn);
        p = Anufft.p;
        Kd = Anufft.Kd;
        Nd = Anufft.Nd;
        Asense = single(cat(4,A.sensmaps{:}));
        
        parfor k=1:2%A.numCoils
            %Q = Q + nufft_adj(B((k-1)*tra_length+1:k*tra_length), A.nufftStruct) .* conj(A.sensmaps{k});
            %Q = Q + nufft_adj(Btemp(:,k), Anufft) .* conj(Asense(:,:,:,k));
            nufft_adj_local(Btemp(:,k), sn, p, Kd, Nd);
        end
        Q = Q / sqrt(prod(A.imageDim));

    else
        Qtmp = zeros(tra_length, A.numCoils);
        Q = zeros(tra_length*A.numCoils, 1);
        parfor k=1:A.numCoils
            %Q((k-1)*tra_length+1:k*tra_length) = nufft((B.*A.sensmaps{k}), A.nufftStruct);
            Qtmp(:,k) = nufft((B.*A.sensmaps{k}), A.nufftStruct);
        end
        Q = reshape(Qtmp,[tra_length*A.numCoils, 1]);
        Q = Q / sqrt(prod(A.imageDim));
        
    end
    
% now B is the operator and A is the vector
else
    Q = mtimes(B',A')';
    
end


function x = nufft_adj_local(y, scaling_factor, interpolation_matrix, Kd, Nd)

% usage:
%   y = nufft_adj_simple(x, scaling_factor, interpolation_matrix, Kd, Nd)
%
% This function returns the same result as nufft_adj.m (although not
% all features are supported), the difference is that the necessary
% data to calculate the result can be passed independently.  The
% arguments can be gained by calling nufft_init.m:
%
% if  nstr = nufft_init(om, Nd, Jd, Kd, nufft_shift, ...), then
% scaling_factor = nstr.sn
% interpolation_matrix = nstr.p
% Kd = nstr.Kd
% Nd = nstr.Nd

dims = size(y);

x = full(interpolation_matrix' * y);
x = reshape(x, [Kd 1]);
x = prod(Kd) * col(ifftn_fast(x));
x = reshape(x, [Kd 1]);

% eliminate zero padding from ends
if length(Nd) == 1
	x = x(1:Nd(1),:);
elseif length(Nd) == 2
	x = x(1:Nd(1),1:Nd(2),:);
elseif length(Nd) == 3
	x = x(1:Nd(1),1:Nd(2),1:Nd(3),:);
else
	error 'only up to 3D implemented currently'
end

x = reshape(x, [prod(Nd) 1]);
x = x .* conj(col(scaling_factor));
x = reshape(x, [Nd dims(2:end)]);
