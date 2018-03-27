function  A = A_opt(sense_map, image_dim, coil_nr, shots_nr, mask)

% usage:
% 
% 
% mask in size [ny nz ns]
% sensemaps in size of [ny nz nc]
% image_dim =  [ny nz]


if nargin==0 % default constructor
    a.adjoint = 0;
    a.sens = [];
    a.image_dim = [];
    a.numCoils = [];
    a.numShots = [];
    a.mask = [];
    a.kspa_dim = [];

else
    a.adjoint = 0;
    a.sens = sense_map;
    a.image_dim = image_dim;
    a.numCoils = coil_nr;
    a.numShots = shots_nr;
    a.mask = mask;
    a.kspa_dim = [size(mask,1) size(mask,2)];
     
end

A = class(a,'A_opt');
