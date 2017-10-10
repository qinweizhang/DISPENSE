function [sense_map, varargout] = calc_sense_map(data_fn, sense_ref_fn, coil_survey_fn, recon_dim,methods)
% calc sens_maps


%get k space

if(strcmp(methods, 'ecalib'))
    
    sense_map = get_sense_map_ecalib(data_fn,recon_dim);
    
elseif (strcmp(methods, 'external'))
    
    [sense_map, sense_Psi] = get_sense_map_external(sense_ref_fn, data_fn, coil_survey_fn, recon_dim);
end

%normalize sense maps
sense_map = normalize_sense_map(sense_map);

if(nargout == 2)
    varargout{1} = sense_Psi;
end


end
