%calculate input phase error

 clear; clc; close all
 load('error_simulation_input.mat');
 
 dummy_shot_nr = 0; % dummy shots used before acquisition starts
 shot_nr_total = 500; 
 phase_error_input_resized = repmat(phase_error_input, [6, 1]); %make sure the length is enough
 ref_shot = 4; % reference shot used in the measured navigator analysis (this shot have no phase error)
 %============================================
 %             global phase error
 %============================================

 phase_step = 0;
 global_phase_error_input = phase_error_input_resized(:,4)/180*pi * phase_step; %in rad 

 %-------compare with navigator fitted global error
 global_phase_error_input_match_rad = global_phase_error_input(dummy_shot_nr+1: (dummy_shot_nr + shot_nr_total)); %remove first 2 dummy shots
 
 i = sqrt(-1);
 measured_global_phase_error = angle(exp(i * global_phase'/2));  % navigator measured phase error has an effect of 2 and maybe 2pi jumpping
 
 figure(1)
 plot((global_phase_error_input_match_rad-global_phase_error_input_match_rad(ref_shot)),'r','LineWidth',2)%the 4th shot is the reference
 hold on; 
 plot(measured_global_phase_error,'b');  
 title('0^t^h order phase error'); xlabel('shots'); ylabel('rad')
 legend('input', 'measured (all channels)')
 
 
 %============================================
 %             linear phase error
 %============================================
 
 phase_step_x = 0.001; GR_x_teff = 0.5172; %in ms
 phase_step_y = 0.001; GR_y_teff = 0.5172; %in ms
 phase_step_z = 0.001; GR_z_teff = 0.5172; %in ms
 
 GR_m = 16.2729; AQ_interval = 0.0066;
 delta_kx_aera_os = GR_m*AQ_interval; % TSE GRm aera change for two consecutive kx points in(ms*mT/m)
 delta_kx_aera = delta_kx_aera_os * 2; % TSE GRm aera change for two consecutive kx points in(ms*mT/m) (remove the oversampled point)
 GR_py_stp = 0.4128; GR_py_teff = 0.5172;
 delta_ky_aera =  GR_py_stp*GR_py_teff; % TSE GRpy aera change for two consecutive ky points in(ms*mT/m)
 GR_pz_stp = 0.4128; GR_pz_teff = 0.5172;
 delta_kz_aera =  GR_pz_stp*GR_pz_teff; % TSE GRpy aera change for two consecutive kz points in(ms*mT/m)

   
 linear_phase_error_input_gr_factor = phase_error_input_resized(dummy_shot_nr+1: (dummy_shot_nr + shot_nr_total),1:3); %remove the first two dummy shots
 
 %Calculate the gradient aera; which is porpotion to kspace displacement in (1/m)
 linear_phase_input_error_x = linear_phase_error_input_gr_factor(:,1) * phase_step_x * GR_x_teff; 
 linear_phase_input_error_y = linear_phase_error_input_gr_factor(:,2) * phase_step_y * GR_y_teff; 
 linear_phase_input_error_z = linear_phase_error_input_gr_factor(:,3) * phase_step_z * GR_z_teff; 
 
 %Calculate the trajecotry offset effect caused by lienar phase error in unit of (k space pixel)
 linear_phase_error_input_x_in_TSEksp_pixel = linear_phase_input_error_x / delta_kx_aera_os;   
%  linear_phase_error_input_x_in_TSEksp_pixel = linear_phase_input_error_x / delta_kx_aera;   
 linear_phase_error_input_y_in_TSEksp_pixel = linear_phase_input_error_y / delta_ky_aera;
 linear_phase_error_input_z_in_TSEksp_pixel = linear_phase_input_error_z / delta_kz_aera;
 
  %-------compare with navigator fitted global error
  
  %in unit of (k space pixel): assumption FOVnav = FOVima; otherwise
  %measured_linear_phase_error_x = measured_linear_phase_error_x * (FOVima / FOVnav)
  nav_kx_dim = 24;    nav_ky_dim = 24;
  measured_linear_phase_error_x = squeeze(linear_phase_xy(1,:,:))* nav_kx_dim / 2 / pi; 
  measured_linear_phase_error_y = squeeze(linear_phase_xy(2,:,:))* nav_ky_dim / 2 / pi; 

  figure(2);
  plot(linear_phase_error_input_x_in_TSEksp_pixel - linear_phase_error_input_x_in_TSEksp_pixel(ref_shot),'r','LineWidth',2)
  hold on
  plot(measured_linear_phase_error_x', 'b');
  title('1^s^t order phase error in x'); xlabel('shots'); ylabel('offset in kspa pixels')
  legend('input', 'measured (all channels)') 
  
  figure(3);
  plot(linear_phase_error_input_y_in_TSEksp_pixel - linear_phase_error_input_y_in_TSEksp_pixel(ref_shot),'r','LineWidth',2)
  hold on
  plot( measured_linear_phase_error_y', 'b');
  title('1^s^t order phase error in y'); xlabel('shots'); ylabel('offset in kspa pixels')
  legend('input', 'measured (all channels)')
  
  save('input_phase_error.mat','linear_phase_error_input_x_in_TSEksp_pixel','linear_phase_error_input_y_in_TSEksp_pixel',...
      'linear_phase_error_input_z_in_TSEksp_pixel', 'global_phase_error_input_match_rad','-append');
  