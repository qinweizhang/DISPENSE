function Q = mtimes(A,B)

% nufft_simple(B .* A.sensmaps{k}, A.scaling_factor, A.interpolation_matrix{m}, A.oversampling)
% nufft_adj_simple(B .* A.sensmaps{k}, A.scaling_factor, A.interpolation_matrix{m}, A.oversampling, A.imageDim)


if strcmp(class(A),'orc_segm_nuFTOperator')
    
    if A.adjoint
        Q = zeros(size(A));
        for k=1:A.numCoils
            tmp = B((k-1)*A.trajectory_length+1:k*A.trajectory_length);

            y = zeros(size(A));
            for m=1:A.nSegments+2
                x = nufft_adj_simple(A.segment_filter{m} .* tmp(A.segment_index{m}), A.scaling_factor, A.interpolation_matrix{m}, A.oversampling, A.imageDim);
                y = y + exp(-1i*A.phasemap*A.t(m)) .* x;
            end
            Q = Q + (y .* conj(A.sensmaps{k}));
        end
        Q = Q / sqrt(prod(A.imageDim));
        

    else
        Q = zeros(A.trajectory_length*A.numCoils, 1);
        for k=1:A.numCoils
            tmp = B .* A.sensmaps{k};

            y = zeros(A.trajectory_length,1);
            for m=1:A.nSegments+2
                x = exp(1i*A.phasemap*A.t(m)) .* tmp;
                y(A.segment_index{m}) = y(A.segment_index{m}) + A.segment_filter{m} .* nufft_simple(x, A.scaling_factor, A.interpolation_matrix{m}, A.oversampling);
            end
            
            Q((k-1)*A.trajectory_length+1:k*A.trajectory_length) = y;
        end
        Q = Q / sqrt(prod(A.imageDim));
        
    end
    
    
% now B is the operator and A is the vector
else
    Q = mtimes(B',A')';
    
end