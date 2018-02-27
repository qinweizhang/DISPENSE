function  A = tSegNUFTOperator(b0_maps, time_segments_par, samples_idx, nuFTop, image_dim, ch_nr)

% usage:
% A = tSegNUFTOperator(b0_maps_resized, nuFTop)
%
% b0_maps_resized:      in size of[x y z]     B0 offresonance map
% time_segments_par:    parameter structure
% samples_idx:          selected recon time point index
% nuFTop:               a defualt nuFFT operator defined by nuFTOperator.m
% image_dim:            recon image dimension
% ch_nr:                channel nr.
%
% IMPLEMENTATION of "Fast, Iterative Imaging Reconstruction for MRI in the Presence of Field Inhomogeneities" Bradley Sutton et al. 
%
% (C) q.zhang Amsterdam, 2017


if nargin==0 % default constructor
    
    s.b0_maps_mHz = [];
    s.tSeg = [];
    s.samples_idx = [];
    s.interpolator = [];
    s.nuFTop = [];
    s.adjoint = [];
    s.ch_nr = [];
    s.image_dim = [];
  
else
    
    s.adjoint = 0;
    s.b0_maps_mHz = b0_maps / 1000;  % % convert Hz(1/s) to mHz(1/ms)
    s.tSeg = time_segments_par;
    s.samples_idx = samples_idx;
    s.nuFTop = nuFTop;
    [s.interpolator, s.tSeg] = tSeg_nuuft_interpolator_init(s.tSeg, s.samples_idx, s.b0_maps_mHz);
    s.ch_nr = ch_nr;
    s.image_dim = image_dim;

end

A = class(s,'tSegNUFTOperator');
