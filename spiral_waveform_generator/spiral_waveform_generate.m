function spiral_waveform_generate
%==========const.=================
max_sample = 18000; %Set in the PPE
gamma = 42577.46778; %Hz/mT
gr_dwell = 0.0064; %in ms
gr_max_str = 40; % in mT/m
gr_max_slewRate = 80; %in mT/m/s

%=========Scan Parameters============
Par.FOV_xy = 0.210;%FOV in M/P in m
Par.FOV_z = 0.040;  %FOV in S in m
Par.VOX_xy = 0.010; %Voxel size in M/P in m
Par.VOX_z = 0.010;  %Voxel size in S in m
%============END=====================

kz_max = 2 * pi / Par.VOX_z;
kxy_max = 2 * pi / Par.VOX_xy;

fun_Rz = @(kz, kz_max)(1+2 * abs(kz) ./ kz_max); %undersample profile in z (User define)
fun_Rr = @(r, kzy_max)(1.5 + abs(r) ./ kxy_max);  %undersample profile in z (User define)
%% Kz and Rz(undersample profile in z)
delta_kz = 2 * pi / Par.FOV_z;
kz = cat(2, fliplr(-delta_kz:-delta_kz:-kz_max),0:delta_kz:kz_max);
Rz = fun_Rz(kz,kz_max);
%%
%Based on equation 1 - 4 in paper "Single shot whole brain imaging using spherical stack of spirals trajectories"
%Mode: Kxy = r(t) * exp(i * theta(t));
current_t = 0;
current_r = 0;
current_rp = 10;
current_rpp = 0;%gr_max_slewRate;

% while(current_t < 100)
tspan = [current_t:gr_dwell: 100];
y0 = current_r;
%Conditional #1: always use the max GRxy strength; expression is explict
[t, y]  = ode45(@ODE_gr_str_condition,tspan,y0);
figure(1); plot(t, y);
%Conditional #2: always use the max GR slew rate; expression is implict

yy0 = [current_r; current_rp];
yyp0 = [current_rp; current_rpp];
[tt, yy]  = ode15i(@ODE_gr_slewrate_conditon, tspan, yy0, yyp0);
figure(2); hold on; plot(tt,yy);






%Conditional #2: always use the max GR strength
%Gxy_slewrate = exp(i * theta) / gamma (i *(r''^2)*(FOVxy/Rr)*(2+i*r*FOVxy/Rr-(r/Rr)*(dRr/dr))+r''*(1+i*r(FOVxy/Rr)));
% i = sqrt(-1);
% kxy_max = 2 * pi / Par.VOX_xy;
% fun_Gxy_str = @(r_1st_dev, r, Rr)(abs(1 / gamma * r_1st_dev * (1 + i * r * Par.FOV_xy / Rr)));
% fun_Gxy_slewrate = @(r_2nd_dev, r_1st_dev, r, Rr, dRr_dr)(abs(1 / gamma*(i * (r_1st_dev.^2)*(Par.FOV_xy/Rr)*(2+i*r*Par.FOV_xy/Rr-(r/Rr)*dRr_dr)+r_2nd_dev*(1+i*r(Par.FOV_xy/Rr)))));
%
% %calculation for r
% stop = false;
%
% %initial values
% r(1) = 0;
% Rr(1) = fun_Rr(abs(r(1)+i*kz(5)), kxy_max);
% r_1st_dev(1) = gr_max_str ./ abs(1 + i * r(1) * Par.FOV_xy/Rr) * gamma;
% r_2nd_dev(1) = gr_max_slewRate * gamma -



    function dydt = ODE_gr_str_condition(t,y)
        %Gxy_str = exp(i * theta) / gamma * r' * (1 + i * r * FOVxy / Rr);
        
        Rr = fun_Rr(abs(y+i*kz(5)), kxy_max);
        dydt = gamma * gr_max_str / abs(1 + i * y * Par.FOV_xy / Rr );
    end

    function f = ODE_gr_slewrate_conditon(t, y, yp)
        
        %Gxy_slewrate = abs(exp(i * theta) / gamma (i *(r''^2)*(FOVxy   /Rr)*(2+i*r   *FOVxy/Rr-(r/Rr)*(dRr/dr))+r''  *(1+i*r   *(FOVxy    /Rr))));
        
        %y(1) = r;  yp(1) = r';
        %y(2) = r'; yp(2) = r'';
        Rr = fun_Rr(abs(y(1)+i*kz(5)), kxy_max);
%         Rr_last = fun_Rr(abs(current_r+i*kz(5)), kxy_max);
        %         dRr_dr = (Rr - Rr_last)/(y(1)-current_r);
        dRr_dr = 1.5;
        
        f = [gr_max_slewRate - abs((i .*yp(1).^2.*(Par.FOV_xy/Rr)*(2+i*y(1)*Par.FOV_xy/Rr-(y(1)/Rr)*(dRr_dr))+yp(2)*(1+i*y(1)*(Par.FOV_xy/Rr))))./gamma;...
%             10*yp(2) + 4*y(1) + Rr;
            yp(1) - y(2)];
    end

    function f = my_ODE15i(t, q, dq)
        f = [dq(1)-q(2);10*dq(2) + 4*q(1)];
    end


% dlmwrite('external_spiral_GRwaveforms.dat',waveform,'delimiter','\t');
%
%
% Ta = 70; %ms
% t = [0:gr_dwell:Ta]; %ms
%
% theta = pi .* t ./ sqrt(1 + t ./ Ta);
% radius = rex

end