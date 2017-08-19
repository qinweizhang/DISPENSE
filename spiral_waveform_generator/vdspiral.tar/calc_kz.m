function kz_all = calc_kz(kz_max, FOV_z, Fz_coeff);


kz = 0;
kz_all = 0;
while kz < kz_max
    
    %current FOVz at kz
    FOVz = 0;
    for coeff = 1:length(Fz_coeff)
        FOVz = FOVz+(Fz_coeff(coeff) .* FOV_z)*(kz/kz_max)^(coeff-1);
    end
    
    delta_kz = 1 / FOVz;
    kz = kz + delta_kz;
    kz_all = [kz_all kz];
end

if(mod(length(kz_all),2) == 0) %kz must be at odd number
    kz_all(end) = [];
end

kz_all = cat(2, -1 * fliplr(kz_all(2:end)), kz_all);

end
