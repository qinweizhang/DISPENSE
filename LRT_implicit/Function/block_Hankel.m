function F_m = block_Hankel(m, pars_Hankel)
% 
% Generate the block-Hankel matrix for the input k-space m
% m has the size of [nx, ny, ns]
% based on Mani et al. Multi-shot Sensitivity-Encoded Diffusion Data Recovery Using Structured Low-Rank Matrix Completion (MUSSELS)
% 
% 2018 Q Zhang - AMC Amsterdam

kernel_size = pars_Hankel.kernel_size;
[nx, ny, ~] = size(m);


if(~pars_Hankel.smooth_constrain) %don't use additional smoothness constrain
    H_m = [];
    for x_idx = 1:(nx - kernel_size(1))
        for y_idx = 1:(ny - kernel_size(2))
            kernel_cube = m(x_idx:(x_idx+kernel_size(1)), y_idx:(y_idx+kernel_size(2)), :);
            H_m = cat(1, H_m, kernel_cube(:).');        
        end
    end
    D_m = H_m;
    
    
    
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
    for x_idx = 1:(nx - kernel_size(1))
        for y_idx = 1:(ny - kernel_size(2))
            kernel_cube = dx_m(x_idx:(x_idx+kernel_size(1)), y_idx:(y_idx+kernel_size(2)), :);
            H_dx_m = cat(1, H_dx_m, kernel_cube(:).');        
        end
    end
    H_dy_m = [];
    for x_idx = 1:(nx - kernel_size(1))
        for y_idx = 1:(ny - kernel_size(2))
            kernel_cube = dy_m(x_idx:(x_idx+kernel_size(1)), y_idx:(y_idx+kernel_size(2)), :);
            H_dy_m = cat(1, H_dy_m, kernel_cube(:).');        
        end
    end
    D_m = cat(1, H_dx_m, H_dy_m);
end
