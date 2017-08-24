function Q = mtimes(A,B)

if strcmp(class(A),'nuFTOperator')
    if A.adjoint
        Q = zeros(size(A));
        for k=1:A.numCoils
            Q = Q + nufft_adj(B((k-1)*A.trajectory_length+1:k*A.trajectory_length), A.nufftStruct) .* conj(A.sensmaps{k});
       end
        Q = Q / sqrt(prod(A.imageDim));

    else
        Q = zeros(A.trajectory_length*A.numCoils, 1);
        for k=1:A.numCoils
            Q((k-1)*A.trajectory_length+1:k*A.trajectory_length) = nufft((B.*A.sensmaps{k}), A.nufftStruct);
        end
        Q = Q / sqrt(prod(A.imageDim));
        
    end
    
% now B is the operator and A is the vector
else
    Q = mtimes(B',A')';
    
end