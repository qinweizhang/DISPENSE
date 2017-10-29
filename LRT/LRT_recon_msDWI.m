function P_recon=LRT_recon_msDWI(kspace,sens,params)
% 4D Low-Rank Tensor reconstruction 

% recon is run as: P_recon=LRT_recon(kspace,sens,params)
% params is a struct- should be made with params_init

% inputs:
% kspace: undersampled kspace with dimensions (x,y,coils,param1,param2)
% sens  : sens maps - size (x,y,coils) 
% params: params is a struct- should be made with params_init

% based on IEEE-TMI paper by Jingfei He  ,2016 

% 2017 J Schoormans and Q Zhang - AMC Amsterdam

fprintf('------------------------------ \n')
fprintf('Low-Rank Tensor Reconstruction \n')
fprintf('------------------------------ \n')
fprintf('   J Schoormans & Q Zhang  \n')
fprintf('       Amsterdam 2017      \n \n')


%%
assert(params.Lg<=params.L3*params.L4,'reduce spatial rank!'); 
assert(~xor(isempty(params.nav_estimate_1),isempty(params.nav_estimate_2)),'Input either both or no subspaces!')
if(params.mix_trajectory)  %mix trajectory; kspace is a cell
    res1  = size(sens,1);  res2 = size(sens,2);
    imagesize = [res1, res2, 1, size(kspace,1), size(kspace,2)];
    tensorsize = imagesize; tensorsize(3) = size(sens,3);
    unfoldedIsize=[imagesize(1)*imagesize(2),imagesize(4)*imagesize(5)];                %coil combined
    unfoldedKsize=[];     %coils separate
        
else %normal case; kspace is a 5-D array
    res1=size(kspace,1);
    res2=size(kspace,2);
    tensorsize=size(kspace);
    imagesize=tensorsize; imagesize(3)=1;
    unfoldedIsize=[size(kspace,1)*size(kspace,2),size(kspace,4)*size(kspace,5)];                %coil combined
    unfoldedKsize=[size(kspace,1)*size(kspace,2)*size(kspace,3),size(kspace,4)*size(kspace,5)];     %coils separate
end
%%
% >>>>>>>>>>>>>>>>>>>>RECON FROM HERE<<<<<<<<<<<<<<<<<<<<<<<<<<<

% FIND MASK
if(iscell(kspace))
    mask = [];
    %======Estimate on eigenvectors for SoSNav multi-shot Diffusion====================%
    %navigator dim; first column
    all_sosnav_k_data = cat(3, kspace{:,1});    %size: kx, ch, shots
    all_sosnav_k_data = permute(all_sosnav_k_data, [3 1 2]);%size: shots, kx, ch 
    S = all_sosnav_k_data(:,:);   %size: shots, kx*ch
    % calculate singular value decomposition
    [left_1,eigen_1,~]=svd(S,'econ');
    params.eigenvals_1=diag(eigen_1);  %output for evaluation
    %navigator estimate are first L left singular vectors
    params.nav_estimate_1=left_1(:,1:params.L3);
    
    %b0 dim: first row
    %convert 2D cartesian kspace b0 image to 1D spiral signal as navigator
    im_b0 = ifft2d( kspace{1,2});
    for c =1:size(im_b0, 3)
        spiral_k_b0(:,c) = params.NUFFT_nav_1ch * im_b0(:,:,c);
    end
    nav_k_b0 = kspace{1,1};
    S = cat(1, spiral_k_b0(:).', nav_k_b0(:).');   %size: shots, kx*ch
    % calculate singular value decomposition
    [left_2,eigen_2,~]=svd(S,'econ');
    params.eigenvals_2=diag(eigen_2);  %output for evaluation
    %navigator estimate are first L left singular vectors
    params.nav_estimate_2=left_2(:,1:params.L4);
    %=============================end=============================================%
else %normal
    mask=squeeze(kspace(:,:,1,:,:))~=0;
end

if isempty(params.nav_estimate_1)                  % subspace estimation 
    [Kx1,Ky1,Kx2,Ky2]=findSharedKpoints(mask,params);
    fprintf('Estimating subspaces...\n')
    fprintf('# of shared ky-kz points of dim 1: %i \n',numel(Kx1))
    fprintf('# of shared ky-kz points of dim 2: %i \n',numel(Kx2))
    
    % 2: estimate subspaces: generalized for non-square shared k-points
    nav_parameter_dim1=[];
    for iter=1:length(Kx1)
        nav_parameter_dim1=cat(1,nav_parameter_dim1,(kspace(Kx1(iter),Ky1(iter),:,:,params.subspacedim1)));
    end
    [nav_estimate_1,params.eigenvals_1]= subspace_estimator_multicoil(squeeze(nav_parameter_dim1),params.L3);
    
    nav_parameter_dim2=[];
    for iter=1:length(Kx2)
        nav_parameter_dim2=cat(1,nav_parameter_dim2,(kspace(Kx2(iter),Ky2(iter),:,params.subspacedim2,:)));
    end
    [nav_estimate_2,params.eigenvals_2]= subspace_estimator_multicoil(squeeze(nav_parameter_dim2),params.L4);
    
    
    
    params.nav_estimate_1=nav_estimate_1;
    params.nav_estimate_2=nav_estimate_2;
    params.params.eigenvals_1=params.eigenvals_1;
    params.params.eigenvals_2=params.eigenvals_2;
else
        fprintf('Subspaces user-defined...\n')
end

% 3: initialize other operators
if strcmp(params.sparsity_transform,'wavelet')
    Psi=opWavelet2(res1,res2,'Daubechies'); %wavelet operator (uses SPOT toolbox (+ other dependencies maybe?)
    operatorsize=[96,96]; % not sure - to do
elseif strcmp(params.sparsity_transform,'TV')
    Psi1=opConvolve(res1,res2,[-1 1],[0 0],'cyclic');
    Psi2=opConvolve(res1,res2,[-1 1]',[0 0],'cyclic');
    Psi=[Psi1;Psi2];
    operatorsize=[res1,res2*2];
elseif strcmp(params.sparsity_transform,'TVOP')
    Psi=TVOP;
    operatorsize=[res1,res2*2]; % not sure - to do 
elseif strcmp(params.sparsity_transform,'I')
    Psi=opEye(res1*res2);
    operatorsize=[res1,res2]; % not sure - to do 

else
    error('sparsity_transform not recognized')
end

if params.normalize_sense %find out how it should be done...
    sensvec=reshape(sens,[size(sens,1)*size(sens,2),size(sens,3)]);
    sens1_mag = reshape(vecnorm(reshape(sensvec, [], size(sens,3)).'), [size(sens,1),size(sens,2)]);
    sens_normalized = bsxfun(@rdivide, sens, sens1_mag);
    sens_normalized(find(isnan(sens_normalized))) = 0;
else sens_normalized=sens;
end
F=MCFopClass;
set_MCFop_Params(F,(sens_normalized),[res1,res2],[tensorsize(4),tensorsize(5)], params.NUFFT_nav_sense, params.trj_length);

%4 zero-filled recon

if(iscell(kspace))
    kspace_1 = LRT_cell_2_mat(kspace); %convert to same size and matrxi unfold
    unfoldedKsize = size(kspace_1); clear kspace;
    
    P1_0=F'*(kspace_1);
    P0=reshape(P1_0,imagesize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
    if params.scaleksp
        [kspace,scaling]= scaleksp(kspace,P0); % scale kspace to ensure consistency over params;
        params.Imref=params.Imref./scaling; %scale ref image with same scaling;
        P1_0=F'*(kspace);
        P0=reshape(P1_0,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
    end
    kspace=reshape(kspace_1,tensorsize);

else
    
    P0=F'*(kspace);
    P1_0=reshape(P0,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
    
    if params.scaleksp
        [kspace,scaling]= scaleksp(kspace,P0); % scale kspace to ensure consistency over params;
        params.Imref=params.Imref./scaling; %scale ref image with same scaling;
        P0=F'*(kspace);
        P1_0=reshape(P0,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
    end
    
    kspace_1=reshape(kspace,unfoldedKsize);
end
figure(21); subplot(211); imshow(abs(P0(:,:,1,params.subspacedim1,params.subspacedim2)),[]); axis off; title('zero filled recon of one frame')
figure(21); subplot(212);  imshow(angle(P0(:,:,1,params.subspacedim1,params.subspacedim2)),[]); axis off; title('phase of zero filled recon of one frame')
figure(22);subplot(311); immontage4D(mask,[0 1]); xlabel('Parameter 1'); ylabel('Parameter 2');
figure(22); subplot(312);  immontage4D(squeeze(abs(P0)));
figure(22); subplot(313); immontage4D(squeeze(angle(P0)),[-pi pi]);
set(0,'DefaultAxesColorOrder',jet(max([size(params.nav_estimate_1,2), size(params.nav_estimate_2,2)]))); 
figure(23); subplot(221); plot(abs(params.nav_estimate_1)); colorbar
figure(23); subplot(222); plot(abs(params.nav_estimate_2)); colorbar
figure(23); subplot(212); hold on; plot(params.eigenvals_1./max(params.eigenvals_1(:)),'r'); plot(params.L3,params.eigenvals_1(params.L3)./max(params.eigenvals_1(:)),'ro');...
    plot(params.eigenvals_2./max(params.eigenvals_2(:)),'b');plot(params.L4,params.eigenvals_2(params.L4)./max(params.eigenvals_2(:)),'bo') ;hold off;
title('eigenvalues for the two subspaces (1=red,2-blue)');
drawnow;
%% ALGO 
alpha=params.alpha;
beta=params.beta; 
mu=params.mu; 
lambda=params.lambda; 

%initialize matrices
if params.inspectLg
params.Lg=params.L3*params.L4;end

[Phi,G,C,A,B,Y,Z]= init_G0(P1_0,Psi,params.nav_estimate_1,params.nav_estimate_2,params.Lg);  

if params.inspectLg;
    C_energies=sum(abs(C).^2,2);
    C_energies=C_energies./max(C_energies(:));
    figure(11); plot(C_energies,'ko-');title('rank relative energies for C');
    fprintf('relE %f: \n',C_energies)
    params.Lg=input('Choose spatial rank: ');
    [Phi,G,C,A,B,Y,Z]= init_G0(P1_0,Psi,params.nav_estimate_1,params.nav_estimate_2,params.Lg);  
end


MSE=[]; 
for iter=1:params.niter
    params.iter=iter; 
    fprintf('\n Outer iteration %i of %i \n',params.iter,params.niter)
    MSE=visualize_convergence(params.iter,MSE,G,C,Phi,params.Imref,imagesize,params.x,params.y);
    
    [Ak,lambda]=soft_thresh_A(G,Y,alpha,lambda,Psi,operatorsize,params);                     %15
    [Bk,mu]=soft_thresh_B(C,Z,mu,beta,params);                              %16
    Gk=precon_conj_grad_G(G,C,Ak,Y,alpha,Psi,kspace_1,Phi,F,params);       %17
    Ck=precon_conj_grad_C(Gk,C,Bk,Z,beta,kspace_1,Phi,F,params);           %18
    Yk=Y+alpha*(Ak-Psi*Gk);
    Zk=Z+beta.*(Bk-Ck);
    
    G=Gk; C=Ck; Y=Yk; Z=Zk; %update iteration

    if params.increase_penalty_parameters
    alpha=alpha*1.5; beta=beta*1.5; end;   
end

P_recon=G*C*Phi;
P_recon=reshape(P_recon,imagesize);

