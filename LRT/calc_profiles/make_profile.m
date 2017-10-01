%calculate scan profile for LRT scan 
%where dimension 1 is the TSE echo number and dimension 2 is T2prep
fprintf('--------------------\n ')
fprintf('Making LRT profile \n')
fprintf('--------------------\n \n')

%%%%%%%%%%%% PARAMETERS TO CHANGE%%%%%%%%%%%%%%%%%%%%%%%
nDim1=40; % TSE dimensions/DTI
nDim2=2; % T2-prep 
ky=109; 
kz=29; 

dim1_bigctr=1; % dimension number of fully sampled center (param dimension 1)
dim2_bigctr=1; % dimension number of fully sampled center (param dimension 1)

bigctrsize=2;
smallctrsize=2;
DTI=1; %1=DTI/T2prep - 0: VFA/T2prep (decides ordering of lines) or VFA/DTI
if(DTI)
ETL = 60;
end
%%%%%%%  CHOOSE ONE OF BOTH OPTIONS
%
nr_points =60; 
undersampling=nr_points./(ky*kz);

% or
% undersampling=0.1;    
% nr_points=ceil(undersampling*ky*kz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%i points per frame, an undersampling factor of %i \n',nr_points,undersampling)
%%%%%%%00

% waiting_time=(521-127)e-3; 
% TR=5.21e-3         %TR in ms;
% TR_shot=nDim1*TR+waiting_time; 
TR_shot=2000e-3;

MC_maxiter=10000; 
visualize=1;
radialflag=1; %radial/linear
linearflag=0; % 0 vertical ordering/ 1 horizontal ordering;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculating some params...
nr_centerpoints=(2*bigctrsize+1)^2; %number of k-points in the center squares; 
assert(nr_points>=nr_centerpoints,'fully sampled centers too big relative to undersampling')
if DTI
    nshots =nDim2*nDim1*nr_points/ ETL;
else
    nshots=nDim2*nr_points; 
end
total_time= TR_shot*nshots; %total time in seconds; 

fprintf('Number of shots: %d, TSE number: %d, total time: %d seconds \n',nshots,nDim1,round(total_time))
assert(dim1_bigctr <= nDim1,'dim1_bigctr should be in range [1, nDim1]')
assert(dim2_bigctr <= nDim2,'dim1_bigctr should be in range [1, nDim1]')

%% add random point to every independent k-space
% performs a Monte Carlo simulation s.t. all have equal # of points
mask=zeros(ky,kz,nDim1,nDim2); 
fprintf('Starting Monte Carlo simulation \n')
for dim1=1:nDim1;
    for dim2=1:nDim2;
        fprintf('Dim 1: %d, Dim 2: %d |',dim1,dim2)
        m=zeros(ky,kz); MC_niter=0;
        if dim1==dim1_bigctr || dim2==dim2_bigctr;
            ctrsize=bigctrsize;
            nr_centerpoints=(2*ctrsize+1)^2; %number of k-points in the center squares; 
        else
            ctrsize=smallctrsize;
            nr_centerpoints=(2*ctrsize+1)^2; %number of k-points in the center squares; 

        end
            while sum(m(:))~=nr_points
                m=rand(size(m))>(1-((nr_points-nr_centerpoints)/(ky*kz)));
                m=addCtr(m,ctrsize);

                MC_niter=MC_niter+1;
                if MC_niter>MC_maxiter
                    error('too many MC iterations - check settings')
                end
                
            end
        fprintf(' Monte Carlo iterations: %d \n',MC_niter)
        mask(:,:,dim1,dim2)=m;
    end
end
clear m
if visualize;   figure(1); clf; imshow(reshape(permute(mask(:,:,1:2,1:2),[1 3 2 4]),[ky*2,kz*2]));
 figure(11); clf; immontage4D(mask,[]);
end

%% order and save

makefilename = @(x) [x,num2str(ky),'_',num2str(kz),'_',num2str(nDim1),'_',num2str(nDim2),...
        '_r',num2str(radialflag),'_l',num2str(linearflag),'_bCtr',num2str(bigctrsize),...
        '(',num2str(dim1_bigctr),',',num2str(dim2_bigctr),')_sCtr',num2str(smallctrsize),'_us',num2str(undersampling),'_nrpoints_',num2str(nr_points)];
    

if ~DTI
    profile_order=profile_ordering(mask,radialflag,linearflag,visualize);
    filename=makefilename('LRT_VFA_T2p_'); 
else
    shot_per_frame = nr_points./ETL;
    assert((round(shot_per_frame)==shot_per_frame),'Profiles in every dynamic must be acquired in integer no. of shots! Change nr_points or ETL!');
    profile_order=profile_ordering_DTI_t2prep(mask,radialflag,linearflag,visualize,shot_per_frame);
    filename=makefilename('LRT_DTI_T2p_'); 
end

 
if ispc()
     cd('L:\basic\divi\Ima\parrec\Jasper\profiles_LRT')
else
     cd(['/home/',getenv('USER'),'/lood_storage/divi/Ima/parrec/Jasper/profiles_LRT']);
end
savemask_LRT(profile_order,filename,visualize)
fprintf('Saved as %s \n',filename)
fprintf('in folder: %s \n',pwd)
fprintf('Finished! \n \n \n ')
