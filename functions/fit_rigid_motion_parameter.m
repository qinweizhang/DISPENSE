
function [rot_vec, trans_cons] = fit_rigid_moion_parameter(non_b0_Ks, ref_Ks)

% function [rot_vec, trans_cons] = fit_rigid_moion_parameter(non_b0_phase, ref_phase, non_b0_Ks, ref_Ks)
% Input 
%     non_b0_Ks: to be corrected navigator data (cpx, K space domain);
% 
%     ref_Ks: reference navigtor data (cpx, K space domain);
%     
% Output
%     rot_vec: vector corresponding to rotaion, apply it to the diffusion K
%     space data, shift the k space cordinate. in unit of (1/m)
% 
%     trans_cons: constant for translation correction, apply it to the K space data, multiply exp(i*trans_cons) to the data
% 
% 
%   (c) Kerrry 2016 AMC based on "Motion-Induced Phase Eror Estimation and Correction in 3D Diffusion Tensor Imaging"



    if(sum((size(non_b0_Ks)-size(ref_Ks)).^2) > 0)
        error('K space matrix are not in the same size...');
    end
    
    [x, y, z] = size(non_b0_Ks);
    
    %-------------initilize esitmation -----------------%
    % argmax|Sref(K)|
    Ref_max_ix = find(abs(ref_Ks(:)) == max(abs(ref_Ks(:))));
    Ref_max_kx = (mod(mod(Ref_max_ix,(x*y)),x));
    Ref_max_ky = (ceil(mod(Ref_max_ix,(x*y))/x));
    Ref_max_kz = (ceil(Ref_max_ix/(x*y)));
    Ref_max = ref_Ks(Ref_max_kx, Ref_max_ky, Ref_max_kz);
    % argmax|Sb(K)|
    non_b0_max_ix = find(abs(non_b0_Ks(:)) == max(abs(non_b0_Ks(:))));
    non_b0_max_kx = (mod(mod(non_b0_max_ix,(x*y)),x));
    non_b0_max_ky = (ceil(mod(non_b0_max_ix,(x*y))/x));
    non_b0_max_kz = (ceil(non_b0_max_ix/(x*y)));
    non_b0_max = non_b0_Ks(non_b0_max_kx, non_b0_max_ky, non_b0_max_kz);
    
    FOVxy = 0.233;  %in x/y directions; in m
    a0(1) = (1/FOVxy)*(non_b0_max_kx - Ref_max_kx);  %voxel diamiter in k space is 1/FOVxy(m-1)
    a0(2) = (1/FOVxy)*(non_b0_max_ky - Ref_max_ky);
    a0(3) = (1/FOVxy)*(non_b0_max_kz - Ref_max_kz);
    a0(4) = angle(non_b0_max) - angle(Ref_max);
    a0 = a0';
    disp(['initial values:', num2str(a0(1)),' ', num2str(a0(2)),' ', num2str(a0(3)),' Const.', num2str(a0(4))]);
    %--------------finished----------------------%
    
    %fft recon
    for i =1: x
        for j = 1:y
            for k =1:z
                phase_mask(i,j,k) = (-1)^(i+j+k);
            end
        end
    end
    non_b0_phase = angle(fftshift(phase_mask.*fftn(non_b0_Ks))); 
    non_b0_mag = abs(fftshift(phase_mask.*fftn(non_b0_Ks))); 
    
    ref_phase  = angle(fftshift(phase_mask.*fftn(ref_Ks))); 
    ref_mag  = abs(fftshift(phase_mask.*fftn(ref_Ks))); 
    
    %unwrap phase?
%     ref_phase = cunwrap(ref_phase,struct('RoundK',true, 'maxblocksize',inf));
%     non_b0_phase = cunwrap(non_b0_phase,struct('RoundK',true, 'maxblocksize',inf));
%     
    
    

    %reshape data to 1D
    nav_ima_intensity_cutoff = 1e6; %purpose: not fit based on noisy background
    non_b0_phase_1d = reshape(non_b0_phase, [x*y*z, 1]);
    b0_phase_1d = reshape(ref_phase, [x*y*z, 1]);   
    r_1d = find(non_b0_mag>nav_ima_intensity_cutoff);
    i=sqrt(-1);
    disp([ ' ',num2str(length(r_1d)/x/y/z*100),'% pixels used for fitting...'])
    
    % 3D correction: fitting function
%     xy_res = FOVxy / x;
%     
%     fun = @(a)abs(exp(i.*non_b0_phase(r_1d)) - exp(i.*(...
%         (ref_phase(r_1d))+...
%         a(1)*(mod(mod(r_1d,(x*y)),x)-floor(x/2))*xy_res+... %xy_resm voxel diamiter in m
%         a(2)*(ceil(mod(r_1d,(x*y))/x)-floor(y/2))*xy_res+...
%         a(3)*(ceil(r_1d/(x*y))-floor(z/2))*xy_res+...     % for z_res not same as xy_res; interpolation to the same?
%         a(4))));
    % 2D correction: fitting function
    xy_res = FOVxy / x;
    
    fun = @(a)abs(exp(i.*non_b0_phase(r_1d)) - exp(i.*(...
        (ref_phase(r_1d))+...
        a(1)*(mod(mod(r_1d,(x*y)),x)-floor(x/2))*xy_res+... %xy_resm voxel diamiter in m
        a(2)*(ceil(mod(r_1d,(x*y))/x)-floor(y/2))*xy_res+...
        a(4))));
     %-------------------------------------%
       
    options.Algorithm = 'trust-region-reflective';
    options.MaxFunEvals = 1000;
    
    a = lsqnonlin(fun,a0,[], [], options); 
    
    rot_vec = a(1:3);
    trans_cons = a(4);
    
    %debug------
    % just use initial value
    rot_vec = a0(1:3);
    trans_cons = a0(4);
    
    phase_def = non_b0_phase-ref_phase;
    
    for i =1: x
        for j = 1:y
            for k =1:z
                calibrated_ref(i, j, k) = ref_phase(i, j, k) + [i-floor(x/2) j-floor(y/2) k-floor(z/2)].*xy_res*rot_vec + trans_cons;
            end
        end
    end
    
    cor_phase_def = non_b0_phase - calibrated_ref;
    
    figure(10); 
    subplot(121);montage(permute(phase_def,[1 2 4 3]),'Displayrange',[-1 1]); colormap jet; 
    colorbar; title(['uncorrected  IniRetvec:',num2str(a0(1)),': ',num2str(a0(2)),': ',num2str(a0(3)),': Trans:',num2str(a0(4))]);
    
    subplot(122); montage(permute(cor_phase_def,[1 2 4 3]),'Displayrange',[-1 1 ]); colormap jet
    colorbar; title(['Corrected def. Retvec:',num2str(rot_vec(1)),': ',num2str(rot_vec(2)),': ',num2str(rot_vec(3)),': Trans:',num2str(trans_cons)]);
     
    drawnow();

end