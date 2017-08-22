
r = [0:0.1:1];
F_coeff = [1 -2 -1 1 6]
for t = 1:length(r)
    Fz_coeff_now(t) = 0;
    for coeff = 1:length(F_coeff)
        Fz_coeff_now(t) = Fz_coeff_now(t)+(F_coeff(coeff))*abs(r(t)).^(coeff-1);
    end
end
figure(5);plot(r, Fz_coeff_now);