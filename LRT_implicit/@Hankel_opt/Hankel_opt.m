function  Hankel = Hankel_opt(image_dim, kernel_size, smooth_flag)

% usage:
% 
% 
% mask in size [ny nz ns]
% sensemaps in size of [ny nz nc]
% image_dim =  [ny nz]


if nargin==0 % default constructor
    a.adjoint = 0;
    a.image_dim = [];
    a.kernel_size = [];
    a.smooth_flag = [];
else
    a.adjoint = 0;
    a.image_dim = image_dim;
    a.kernel_size = kernel_size;
    a.smooth_flag = smooth_flag;
     
end

Hankel = class(a,'Hankel_opt');
