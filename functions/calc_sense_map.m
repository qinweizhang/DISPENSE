function [sense_map, varargout] = calc_sense_map(data_fn, sense_ref_fn, coil_survey_fn, recon_dim,methods, varargin)
% calc sens_maps
if(nargin == 7)
    os = varargin{1};
    other_pars = varargin{2};
elseif(nargin == 6)
    os = varargin{1};
    other_pars = [];
else
    os = [1 1];
    other_pars = [];
end

%get k space

if(strcmp(methods, 'ecalib'))
    
    sense_map = get_sense_map_ecalib(data_fn,recon_dim);
    
elseif (strcmp(methods, 'external'))
    
    [sense_map, sense_Psi] = get_sense_map_external(sense_ref_fn, data_fn, coil_survey_fn, recon_dim, os, other_pars);
end

%normalize sense maps
sense_map = normalize_sense_map(sense_map);

if(nargout == 2)
    varargout{1} = sense_Psi;
end


end
