function k_optimized = smooth_traj_end(k)

gr_dwell    = 6.4e-6;   %in s
gamma       = 4257.6;

assert(abs(angle(k(end)))<eps);

g_y = 1/gamma * cat(2, 0, diff(imag(k))./gr_dwell);

g_x = 1/gamma * cat(2, 0, diff(real(k))./gr_dwell);

%put g_x(end) g_y(end) to res * g_end
blip_sr = 9000; %G/cm/s  90 mT/m/ms
res = 0;
slop = (1-res)*max(abs(g_x(end)/blip_sr), abs(g_y(end)/blip_sr)); %seconds
slop = round(slop/gr_dwell)*gr_dwell;
blip_sr_x = -g_x(end)/slop;
blip_sr_y = -g_y(end)/slop;

str_step_x = blip_sr_x * gr_dwell;
g_x_tail = [g_x(end):str_step_x:res*g_x(end)];
g_x_tail = g_x_tail(2:end);

str_step_y = blip_sr_y * gr_dwell;
g_y_tail = [g_y(end):str_step_y:res*g_y(end)];
g_y_tail = g_y_tail(2:end);

g_x_final = [g_x g_x_tail];
g_y_final = [g_y g_y_tail];


kx = zeros(1, length(g_x_final));
ky = zeros(1, length(g_x_final));
for idx = 2:length(g_x_final)
    kx(idx) = kx(idx-1) + gamma*gr_dwell*g_x_final(idx);
    ky(idx) = ky(idx-1) + gamma*gr_dwell*g_y_final(idx);
end
k_optimized = kx + i*ky;

phi_end = angle(k_optimized(end));
k_optimized =  k_optimized .* exp(-i.*phi_end); %always ends at 0 degree

end