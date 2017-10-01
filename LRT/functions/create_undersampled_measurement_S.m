% making undersampled k-space measurement
% with fully sampled navigator k-space centers 

% TO DO: -add var density and more us options
rng('default') 
rng(3); % random seed to control phantom generation
res=128;


%chose simulation m
switch simulation_typ
    case 'qMATCH'
        %>>>>>>>>>>>>>> generate VFA TSE & T2 weighted phantoms<<<<<<<<<<<<<<<<<<<<
        T1vals=[500 400 300 200 1000 500 400 300 200 1000].*1e-3; % in seconds
        T2vals=[20 30 40 50 50 55 50 60 70 10].*1e-3; % in seconds
        
        T2prep=[10:10:97].*0.001; % in seconds
        TI=1e-3.*[100:100:2000];
        
        I=T2_T1_phantom(res,T2vals,T1vals,T2prep,TI,1,complexsim);
        
        if complexsim
            fprintf('adding 8pi phase gradient over image \n')
            phase_gradient=exp(-1i*linspace(0,8*pi,res));
            I=bsxfun(@times,I,phase_gradient);
        end
    
    
    case 'DWI_n_T2W'
        %>>>>>>>>>>>>>>> generate Diffusion & T1 (T2?) weigthed phantoms<<<<<<<<<<<
        ADCvals=[1:5].*1e-3;
        T1vals=[20 50 100 200 800].*1e-3;
        
        TE=[10:10:100].*0.001; % in seconds
        bvals=[0:100:900]; %units?
        
        I=diffusion_T1_phantom(res,ADCvals,T1vals,TE,bvals,1);
        
    case 'VFA_n_T2W'
        %>>>>>>>>>>>>>> generate VFA TSE & T2 weighted phantoms<<<<<<<<<<<<<<<<<<<<
        T1vals=[500 400 300 200 1000 500 400 300 200 1000].*1e-3; % in seconds
        T2vals=[20 30 40 50 50 55 50 60 70 10].*1e-3; % in seconds
        
        T2prep=[10:10:97].*0.001; % in seconds
        TEes=4.*0.001; % in seconds
        ETL = 20;
        
        I=VFA_TSE_T2_T1_phantom(res,T2vals,T1vals,T2prep,TEes,ETL,1,complexsim);
        
        if complexsim
            fprintf('adding 8pi phase gradient over image \n')
            phase_gradient=exp(-1i*linspace(0,8*pi,res));
            I=bsxfun(@times,I,phase_gradient);
        end
        
    case 'DTI_n_T2W'
        %>>>>>>>>>>>>>> generate DTI & T2 weighted phantoms<<<<<<<<<<<<<<<<<<<<
        T2vals=[20 30 40 50 50 55 50 60 70 10].*1e-3; % in seconds
        %simple DTI model is considered. i.e. no fiber crossing
        %S_bj = S0 .* exp(-b .* Gj * D * Gj');
        %Gj = [x, y, z]; jth diffusion gradient direction in scanner coordinate
        %D =  [Dxx Dxy Dxz; Dyx Dyy Dyz; Dzx Dzy Dzz]
        %Dxx, Dyy, Dzz: diffusion coefficient in x, y, z direction, in scanner coordinate
        %Dxy, Dyz, Dxz: correlation of diffusion in any two direction
        
        Dxx = [1.5 1.0 1.9 1.7 2.5 0.8 1.6 0.6 2.1 2.9].*1e-3; 
        Dyy = [0.3 0.6 0.5 2.1 3.0 1.7 0.5 1.3 0.5 2.9].*1e-3; 
        Dzz = [0.7 0.5 0.7 0.6 0.5 0.7 1.2 1.4 1.6 2.9].*1e-3; 
        Dxy = [1.5 1.0 1.9 1.7 2.5 0.8 1.6 0.6 2.1 2.9].*1e-3; 
        Dxz = [1.5 1.0 1.9 1.7 2.5 0.8 1.6 0.6 2.1 2.9].*1e-3; 
        Dyz = [1.5 1.0 1.9 1.7 2.5 0.8 1.6 0.6 2.1 2.9].*1e-3; 
        D = reshape(cat(1,Dxx, Dxy, Dxz, Dxy, Dyy, Dyz, Dxz, Dyz, Dzz), 3, 3, length(Dxx));
        
        
        T2prep=[10:10:97].*0.001; % in seconds
        b_value = 600; %in s/mm2
        diffusion_gradient = search_DTI_gr_table('IDIFF_SCHEME_OPT12_OVERPLUS');
        
        visualize = 1;
        I=DTI_T2_phantom(res,T2vals,D,b_value,diffusion_gradient,T2prep,visualize,complexsim);
     
    otherwise
        error('Simulation type is unknown...')
end

% give background signal testing sense estimation
% I = I + 0.1;
% end

figure(1); Q=[];
for ii=1:size(I,3)
    J=[];
    for jj=1:size(I,4);
        J=[J,abs(I(:,:,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[0 1])
xlabel('Prameter 1')
ylabel('Prameter 2')
clear Q J 
%%  add coil sensitivity information 
[Ic,sens]=addcoilsensitvity_to_simulated(I,ncoils);

%% to k space 
d=fftshift(fftshift(fft(fft(ifftshift(ifftshift(Ic,2),1),[],1),[],2),1),2)./sqrt(res*res);
figure(2); imshow(abs(d(:,:,1,1)),[]); axis off; title('k-space')

%% make undersampling mask  (independent of coil dimensions!)
params = initprofile_params()
params.undersampling=uf;

params.ky=res; params.kz=res; 
params.nDim1=size(I,3); params.nDim2=size(I,4); 
params.bigctrsize=10;

mask= make_profile_function(params);
mask=permute(mask,[1 2 5 3 4]); %such that dimensions are ok
%% make undersampled measurement; 
du=bsxfun(@times,mask,d); 

%add noise
if noiselevel>0
du=du+(randn(size(du)).*mean(du(:)).*noiselevel).*(du~=0);
end









