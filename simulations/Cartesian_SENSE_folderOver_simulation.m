%% 2D SENSE map 128*128

% generate 128*128 cartesian trajectory
k_cartesian = bart('phantom -s8 -x128 -k');

% use ecalib to get sens maps
bart_com_2 = sprintf('resize -c 0 %d 1 %d',128,128);
k_cartesian_rs = bart(bart_com_2, k_cartesian);

sens_maps_cartesian = bart('ecalib -m1', k_cartesian_rs);
 
%% This gives the expected fold-over image and correct sense recon

% generate 64*64 cartesian trajectory
cart_traj = bart('traj -x 64 -y 64');

% make gaps in the trajectory
K_space_ideal = bart('phantom -s8 -k -t', cart_traj * 2);

% Use the trajectory without gaps 
bart_com_1 = sprintf('nufft -i -l0.01 -d%d:%d:1 ',64,64); %recon size doesn't matter; only for interpolation 
ima_nufft = bart(bart_com_1, bart('scale 1',cart_traj) , K_space_ideal);

figure; imshow(abs(bart('rss 8',ima_nufft)),[])
figure; montage(permute(abs(ima_nufft),[1 2 3 4]),'displayrange',[]); colormap gray

%SENSE recon: use the trajecotry with gaps & sens maps with desired image size
ima_trj_pics = squeeze(bart('pics -S -r0.001 -t', ...
    cart_traj*2, K_space_ideal, sens_maps_cartesian)); %or sens_maps_cartesian
figure; imshow(abs(ima_trj_pics),[])

%% NUFFT simulation
%trajectory should scale to [-pi pi] = bart_nufft interval = 1; 
%scale to e.g. [-2pi 2pi] = bart_nufft interval = 2  >>> big FOV and image repetition
cart_traj_nufft = reshape(cart_traj, 3, 64*64)';
cart_traj_nufft = pi/32*cart_traj_nufft;

%correct direct recon
 A=nuFTOperator(cart_traj_nufft(:,1:2),[64, 64],ones(64, 64,8),6); 
 im_recon_nufft=regularizedReconstruction(A,col(K_space_ideal(1,:,:,:)),@L2Norm,0.5,'maxit',25);
       
 figure; imshow(abs(im_recon_nufft),[])      
 
 %correct sense recon
  A=nuFTOperator(cart_traj_nufft(:,1:2),[128, 128],squeeze(sens_maps_cartesian),6); 
 im_recon_nufft=regularizedReconstruction(A,col(K_space_ideal(1,:,:,:)),@L2Norm,0.5,'maxit',25);
       
 figure; imshow(abs(im_recon_nufft),[])    