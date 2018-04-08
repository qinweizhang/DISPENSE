function res = mtimes(a,b)
%
% Generate the block-Hankel matrix for the input k-space m
% m has the size of [nx, ny, ns]
% based on Mani et al. Multi-shot Sensitivity-Encoded Diffusion Data Recovery Using Structured Low-Rank Matrix Completion (MUSSELS)
%
% 2018 Q Zhang - AMC Amsterdam

kernel_size = a.kernel_size;
nx = a.image_dim(1);
ny = a.image_dim(2);

if a.adjoint  % H_m --> m
    %% H_m --> m
    shots = size(b, 2) / prod(kernel_size); 
    m_mask = zeros(nx, ny);
    m = zeros(nx, ny, shots);
    
    if(~a.smooth_flag) %don't use additional smoothness constrain
        for idx = 1:size(b, 1)
            [x, y] = ind2sub([nx-kernel_size(1)+1, ny-kernel_size(2)+1], idx);
            kernel_cube = reshape(b(idx, :), [kernel_size, shots]);
            m(x:x+kernel_size(1)-1, y:y+kernel_size(2)-1,:) =  m(x:x+kernel_size(1)-1, y:y+kernel_size(2)-1,:) + kernel_cube;
            m_mask(x:x+kernel_size(1)-1, y:y+kernel_size(2)-1) = m_mask(x:x+kernel_size(1)-1, y:y+kernel_size(2)-1) + ones(kernel_size); 
        end
        res = bsxfun(@rdivide, m, m_mask);
    else %use additional smoothness constrain
        %TODO
    end
    
    
    
    
else  %m --> H_m
    %% m --> H_m
    [nx_m, ny_m, ~] = size(b); 
    assert( nx_m == nx);
    assert( ny_m == ny);
    m = b;
    
    if(~a.smooth_flag) %don't use additional smoothness constrain
        H_m = [];
        for y_idx = 1:(ny - kernel_size(2)+1)
            for x_idx = 1:(nx - kernel_size(1)+1)        
                kernel_cube = m(x_idx:(x_idx+kernel_size(1)-1), y_idx:(y_idx+kernel_size(2)-1), :);
                H_m = cat(1, H_m, kernel_cube(:).');
            end
        end
        res = H_m;
        
        
        
    else %use additional smoothness constrain
        i = sqrt(-1);
        dx_m = m;
        for x_idx = 1:nx
            kx = x_idx / nx - 0.5; % How to scale kx? [-0.5, 0.5]??????
            dx_m(x_idx,:) = dx_m(x_idx,:) .* (-i * 2 * pi * kx);
        end
        dy_m = m;
        for y_idx = 1:ny
            ky = y_idx / ny - 0.5; % How to scale ky? [-0.5, 0.5]??????
            dy_m(:, y_idx) = dy_m(:, y_idx) .* (-i * 2 * pi * ky);
        end
        
        
        H_dx_m = [];
        for y_idx = 1:(ny - kernel_size(2)+1)
            for x_idx = 1:(nx - kernel_size(1)+1)
                kernel_cube = dx_m(x_idx:(x_idx+kernel_size(1)-1), y_idx:(y_idx+kernel_size(2)-1), :);
                H_dx_m = cat(1, H_dx_m, kernel_cube(:).');
            end
        end
        H_dy_m = [];
        for y_idx = 1:(ny - kernel_size(2)+1)
            for x_idx = 1:(nx - kernel_size(1)+1)
                kernel_cube = dy_m(x_idx:(x_idx+kernel_size(1)-1), y_idx:(y_idx+kernel_size(2)-1), :);
                H_dy_m = cat(1, H_dy_m, kernel_cube(:).');
            end
        end
        res = cat(1, H_dx_m, H_dy_m);
    end
    
end




