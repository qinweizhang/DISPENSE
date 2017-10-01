function params = initprofile_params()
%make parameter structure for profile making function
% adds initial values - to be changed by user before calling!!


params=struct; 
params.nDim1=10; % TSE dimensions/DTI
params.nDim2=10; % T2-prep 
params.ky=100; 
params.kz=100; 
params.dim1_bigctr=1; % dimension number of fully sampled center (param dimension 1)
params.dim2_bigctr=1; % dimension number of fully sampled center (param dimension 1)
params.bigctrsize=1;
params.smallctrsize=0;
params.DTIflag=0; %1=DTI/T2prep - 0: VFA/T2prep (decides ordering of lines)
params.ETL=69;
params.visualize=1;
params.radialflag=1; %radial/linear
params.linearflag=0; % 0 vertical ordering/ 1 horizontal ordering;
params.undersampling=[]; %
params.vardens=1;


params.TR_shot=1; 
end