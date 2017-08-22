function kz_connect = connect_kz(kz_target, kz_now )

gr_dwell    = 6.4e-6;   %in s
gamma       = 4257.6;

d_kz = kz_target - kz_now;
if(d_kz>=0)
    blip_sr = 2000; %G/cm/s  20 mT/m/ms
else
    blip_sr = -2000; %G/cm/s  20 mT/m/ms  
end
slop = sqrt(d_kz / gamma /blip_sr); %seconds
slop = round(slop/gr_dwell)*gr_dwell;
blip_sr = d_kz/gamma/slop^2;

gz = [0:gr_dwell:slop].*blip_sr;
gz = [gz fliplr(gz(1:end-1))];
kz_connect = ones(1, length(gz)).*kz_now;
for idx = 2:length(gz)
    kz_connect(idx) = kz_connect(idx-1) + gamma*gr_dwell*gz(idx);
end

end