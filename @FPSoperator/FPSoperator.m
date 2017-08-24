function  A = FPSoperator(sensemaps, phase_error_maps, image_dim, coil_nr, shots_nr, mask)

% usage:
% 
% 
% 
% sensemaps in size of [nx ny nz nc]
% phase error maps in size of [nx ny nz nc nshots]
% image_dim =  [nx ny nz]


if nargin==0 % default constructor
    a.adjoint = 0;
    a.sens = [];
    a.pe = [];
    a.image_dim = [];
    a.numCoils = [];
    a.numShots = [];
    a.mask = [];

else
    a.adjoint = 0;
    a.sens = sensemaps;
    a.pe = phase_error_maps;
    a.image_dim = image_dim;
    a.numCoils = coil_nr;
    a.numShots = shots_nr;
    a.mask = mask;
     
end

A = class(a,'FPSoperator');
