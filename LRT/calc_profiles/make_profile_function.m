%calculate scan profile for LRT scan 
%where dimension 1 is the TSE echo number and dimension 2 is T2prep
function mask= make_profile_function(params)

fprintf('--------------------\n ')
fprintf('Making LRT profile \n')
fprintf('--------------------\n \n')
  
%%%%%%%%%%%% PARAMETERS TO CHANGE%%%%%%%%%%%%%%%%%%%%%%%
nDim1=params.nDim1; % TSE dimensions/DTI
nDim2=params.nDim2; % T2-prep 
ky=params.ky; 
kz=params.kz; 
dim1_bigctr=params.dim1_bigctr; % dimension number of fully sampled center (param dimension 1)
dim2_bigctr=params.dim2_bigctr; % dimension number of fully sampled center (param dimension 1)

bigctrsize=params.bigctrsize;
smallctrsize=params.smallctrsize;

DTI=params.DTIflag; %1=DTI/T2prep - 0: VFA/T2prep (decides ordering of lines)
ETL = params.ETL;


%%%%%%%  CHOOSE ONE OF BOTH OPTIONS
if isempty(params.undersampling)
    error('please specify undersampling: either as a fraction (0,1] or as the number of points sampled per frame >1')
elseif params.undersampling>1
    nr_points=params.undersampling;
    undersampling=nr_points./(ky*kz);
else
undersampling=params.undersampling;    nr_points=ceil(undersampling*ky*kz);
end
fprintf('%i points per frame, an undersampling factor of %i \n',nr_points,undersampling)
%%%%%%%00

% waiting_time=(521-127)e-3; 
% TR=5.21e-3         %TR in ms; 
% TR_shot=nDim1*TR+waiting_time; 

MC_maxiter=10000; 

visualize=params.visualize;
radialflag=params.radialflag; %radial/linear
linearflag=params.linearflag; % 0 vertical ordering/ 1 horizontal ordering;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculating some params...
nr_centerpoints=(2*bigctrsize+1)^2; %number of k-points in the center squares; 
assert(nr_points>=nr_centerpoints,'fully sampled centers too big relative to undersampling')
if DTI
    nshots =nDim2*nDim1*nr_points/ ETL;
else
    nshots=nDim2*nr_points; 
end
total_time= params.TR_shot*nshots; %total time in seconds; 

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
if visualize;   
%     figure(1); clf; imshow(reshape(permute(mask(:,:,1:max(2,size(mask,3)),1:max(2,size(mask,4))),[1 3 2 4]),[ky*2,kz*2]));
 figure(11); clf; immontage4D(mask(:,:,1:3,1:3),[]);
end

end
