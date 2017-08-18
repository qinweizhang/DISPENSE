function normed_sense_map = normalize_sense_map(raw_sens_map)

%normallized sense map

assert(length(size(raw_sens_map)) == 4);

sensvec=reshape(raw_sens_map,[size(raw_sens_map,1)*size(raw_sens_map,2)*size(raw_sens_map,3), size(raw_sens_map, 4)]);
sens1_mag = reshape(colnorm2(sensvec.'), [size(raw_sens_map,1),size(raw_sens_map,2), size(raw_sens_map,3)]);
normed_sense_map = bsxfun(@rdivide, raw_sens_map, sens1_mag);

normed_sense_map(find(isnan(normed_sense_map))) = 0;
normed_sense_map(find(isinf(normed_sense_map))) = 0;

end