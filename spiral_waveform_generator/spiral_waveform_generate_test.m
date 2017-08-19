function spiral_waveform_generate_test
clear all; clc; close all;
%==========const.=================
max_sample = 18000; %Set in the PPE
gamma = 1; %42.57746778; %(1/ms/mT)
gr_dwell = 6.4e-3; %in ms
gr_max_str = 40; % in mT/m
gr_max_slewRate = 50; %in mT/m/ms

%=========Scan Parameters============
Par.FOV_xy = 0.210;%FOV in M/P in m
Par.FOV_z = 0.040;  %FOV in S in m
Par.VOX_xy = 0.010; %Voxel size in M/P in m
Par.VOX_z = 0.010;  %Voxel size in S in m
%============END=====================

kz_max = 2 * pi / Par.VOX_z;
kxy_max = 2 * pi / Par.VOX_xy;

fun_Rz = @(kz, kz_max)(1+2 * abs(kz) ./ kz_max); %undersample profile in z (User define)
fun_Rr = @(r, kxy_max)(1.5 + abs(r) ./ kxy_max);  %undersample profile in z (User define)
%% Kz and Rz(undersample profile in z)
delta_kz = 2 * pi / Par.FOV_z;
kz = cat(2, fliplr(-delta_kz:-delta_kz:-kz_max),0:delta_kz:kz_max);
Rz = fun_Rz(kz,kz_max);
%%
%Based on equation 1 - 4 in paper "Single shot whole brain imaging using spherical stack of spirals trajectories"
%Mode: Kxy = r(t) * exp(i * theta(t));
current_t = 0;
current_r = 0;
current_rp = 0;
current_rpp = 0;
current_Gxy = [];
current_Gxy_SR = [];

all_t = [0];
all_r = [0];
all_rp = [0];
all_rpp = [0];
all_Gxy = [];
all_Gxy_SR = [];

stop = false;
updated = false;
while(~stop)
    tspan = [current_t:gr_dwell: 50];
    
    %Conditional #1: always use the max GRxy strength; expression is explict
    y0 = current_r;
    [t r]  = ode45(@ODE_gr_str_condition,tspan,y0);
    
    rp = cat(1, 0, diff(r))./gr_dwell;
    rpp = cat(1, 0, diff(rp))./gr_dwell;
    Gxy = calc_Gxy(r, rp);
    Gxy_SR = calc_Gxy_slewrate(r, rp, rpp);
    
    %Conditional #2: always use the max GR slew rate; expression is implict
    yy0 = [current_r; current_rp];
    yyp0 = [current_rp; current_rpp];
    options = odeset('RelTol',1e1); 
    [tt, rr]  = ode15i(@ODE_gr_slewrate_conditon, tspan, yy0, yyp0, options);
    
%     if(length(rr) == 3)
        rrp = rr(:,2);
        rrpp = cat(1, 0, diff(rrp))./gr_dwell;
        Gxy_2 = calc_Gxy(rr(:,1), rrp);
        Gxy_SR_2 = calc_Gxy_slewrate(rr(:,1), rrp, rrpp);
%     else
%         disp('ode15i for condition #2 failed');
%         rr = zeros(3,2);
%         rrp = zeros(3,1);
%         rrpp = zeros(3,1);
%         Gxy_2 = zeros(3,1);
%         Gxy_SR_2 = zeros(3,1);
%     end
        
%     
%     if(current_t == 0) %first interation. take the first point from codition #2.
%         all_t = [0];
%         all_r = [rr(1)];
%         all_rp = [rrp(1)];
%         all_rpp = [rrpp(1)];
%         all_Gxy = [Gxy_2(1)];
%         all_Gxy_SR = [Gxy_SR_2(1)];
%     end
%     
%     
%     if((Gxy_SR(2)<gr_max_slewRate)||rp(2)<=rrp(2)) %Coditional #1 is more strict
%         all_t = [all_t  current_t + gr_dwell];
%         all_r = [all_r r(2)];
%         all_rp = [all_rp rp(2)];
%         all_rpp = [all_rpp rpp(2)];
%         all_Gxy = [all_Gxy Gxy(2)];
%         all_Gxy_SR = [all_Gxy_SR Gxy_SR(2)];
%         updated = true;
%     else %Coditional #2 is more strict
%         all_t = [all_t  current_t + gr_dwell];
%         all_r = [all_r rr(2)];
%         all_rp = [all_rp rrp(2)];
%         all_rpp = [all_rpp rrpp(2)];
%         all_Gxy = [all_Gxy Gxy_2(2)];
%         all_Gxy_SR = [all_Gxy_SR Gxy_SR_2(2)];
%         updated = true;
%     end
%     
    %update initial value for next iteration
%     current_t = all_t(end);
%     current_r = all_r(end);
%     current_rp = all_rp(end);
%     current_rpp = all_rpp(end);
    
%     if (current_r>kxy_max)
%         disp('r reached kxy_max!');
        stop = true;
%     end
%     
%     if(~updated)
%         disp('both conditions faild!');
%         stop = true;
%     else
%         update = false;
%     end
%     
    if(1)
        figure(1);
        subplot(211); hold on;
        plot(t, r,'g');        plot(t, rp, 'b','LineWidth',2); plot(t, rpp,'r');
%         plot(tt, rr(:,1),'g'); plot(tt, rr(:,2),'b','LineWidth',2); plot(tt, rrpp, 'r');
%         plot(all_t, all_r,'g'); plot(all_t, all_rp,'b','LineWidth',2); plot(all_t, all_rpp, 'r');
        legend('r','rp','rpp','r_2','rp_2','rpp_2');
        
        subplot(212); hold on;
        plot(t, Gxy, 'r');    plot(t, Gxy_SR,'k'); %ylim([0 200])
%         plot(tt, Gxy_2, 'r'); plot(tt, Gxy_SR_2,'k'); %ylim([0 200])
%         plot(all_t, all_Gxy, 'r'); plot(all_t, all_Gxy_SR,'k');
        legend('Gxy', 'Gxy_S_R','Gxy_2','Gxy_S_R_2');
        
        drawnow();
    end
    
    
end
disp('done');



%% local functions

    function dydt = ODE_gr_str_condition(t,y)
        %-----------------------------------------------------------------
        %Gxy_str = exp(i * theta) / gamma * r' * (1 + i * r * FOVxy / Rr);
        %------------------------------------------------------------------
        Rr = fun_Rr(abs(y+i*kz(5)), kxy_max);
        dydt = gamma * gr_max_str ./ abs(1 + i * y * Par.FOV_xy ./ Rr );
    end



    function Gxy = calc_Gxy(r, rp)
        Rr = fun_Rr(abs(r+i*kz(5)), kxy_max);
        Gxy = abs(1/gamma .* rp .* (1 + i .* r .* (Par.FOV_xy./Rr)));
    end




    function f = ODE_gr_slewrate_conditon(t, y, yp)
        %----------------------------------------------------------------
        %Gxy_slewrate = abs(exp(i * theta) / gamma (i *(r'^2)*(FOVxy   /Rr)*(2+i*r   *FOVxy/Rr-(r/Rr)*(dRr/dr))+r''  *(1+i*r   *(FOVxy    /Rr))));
        %y(1) = r;  yp(1) = r';
        %y(2) = r'; yp(2) = r'';
        %-----------------------------------------------------------------
        Rr = 1;%  fun_Rr(abs(y(1)+i*kz(5)), kxy_max);
        %         Rr_last = fun_Rr(abs(current_r+i*kz(5)), kxy_max);
        %         dRr_dr = (Rr - Rr_last)/(y(1)-current_r);
        dRr_dr = 1.5;
        
        f = [%gr_max_slewRate - abs((i .*yp(1).^2.*(Par.FOV_xy/Rr)*(2+i*y(1)*Par.FOV_xy/Rr-(y(1)/Rr)*(dRr_dr))+yp(2)*(1+i*y(1)*(Par.FOV_xy/Rr))))./gamma;...
            gr_max_slewRate * gamma - abs((i .*yp(1).^2.*(Par.FOV_xy/Rr)*(2+i*y(1)*Par.FOV_xy/Rr-(y(1)/Rr)*(dRr_dr))+yp(2)*(1+i*y(1)*(Par.FOV_xy/Rr))));...
            yp(1) - y(2)];
    end




    function Gxy_slewrate = calc_Gxy_slewrate(r, rp, rpp)
        Rr = 1;
        dRr_dr = 1.5;
        Gxy_slewrate = abs(1 / gamma .* (i *(rp.^2)*(Par.FOV_xy  /Rr).*(2+i.*r.*Par.FOV_xy/Rr-(r./Rr).*(dRr_dr))+rpp.*(1+i.*r.*(Par.FOV_xy./Rr))));
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