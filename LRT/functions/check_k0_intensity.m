function check_k0_intensity(MR,par_dim1, par_dim2)

data = MR.Data;

typ = MR.Parameter.Labels.Index.typ;
ky =MR.Parameter.Labels.Index.ky;
kz =MR.Parameter.Labels.Index.kz;
ky_typ1 = ky(find(typ == 1));
kz_typ1 = kz(find(typ == 1));
coil_nr = length(MR.Parameter.Labels.CoilNrs);
ky_typ1_1coil = ky_typ1(1:coil_nr:end);
kz_typ1_1coil = kz_typ1(1:coil_nr:end);

data_coil_nr = MR.Parameter.Recon.ACNrVirtualChannels;
data_1coil = data(:,1:data_coil_nr:end);

k0_locations = find((ky_typ1_1coil == 0).*(kz_typ1_1coil == 0));
k0_data =  data_1coil(:,k0_locations);
k0_data = max(abs(k0_data),[], 1);
figure(210); plot(k0_data);

k0_data = reshape(k0_data, par_dim1, par_dim2);
figure(211); mesh(double(k0_data));
end
