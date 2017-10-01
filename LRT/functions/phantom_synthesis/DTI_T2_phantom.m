%
function I=VFA_TSE_T2_T1_phantom(N,T2vals,D, b_value, diffusion_gradient,T2prep, varargin )
%N : resolution
% T2vals = vector of T2 values
% D = DTI tensor

% diffusion_gradient: vector of GRs
% T2prep

%optional: visualize true/false
%% EXAMPLE 
%I=VFA_TSE_T2_T1_phantom(128,[20 40 60 80 100].*1e-3,[300 500 800 1000 1500].*1e-3,[10:3:97].*0.001,4.*0.001,30,1);

%%
assert(length(T2vals)==size(D,3)); %want same number of bvals as T1vals (= number of phantoms)

if nargin>6
    visualize = varargin{1};
    complexsim=varargin{2};
else
    complexsim=false;
    visualize =false;
end


%% START CODE

for nphantom=1:length(T2vals)
    E(1)=0.4+0.6*rand; %intensity
    E(2)= 0.1+rand*0.2; %length
    E(3)= 0.1+rand*0.2; %width
    E(4)=-0.5+(rand); %x-coord of middle 
    E(5)=-0.5+(rand); %y-coord of middle
    E(6)= rand*360; %angle   
    P{nphantom} = phantom(E,N) ;
    
end




% make simulation images
I=zeros(N,N,length(diffusion_gradient),length(T2prep));
if complexsim
  for nphantom=1:length(T2vals)
    phase(nphantom)=exp(1i*rand*2*pi);
    fprintf('phase phantom %d is %4.2f \n',nphantom,angle(phase(nphantom)))
  end  
end

for ii=1:length(diffusion_gradient)
    for jj=1:length(T2prep)
        
        for nphantom=1:length(T2vals) % for all separate phantoms
            
            if complexsim %add random phase to phantom
                I(:,:,ii,jj)=I(:,:,ii,jj)+(P{nphantom} .* exp(-T2prep(jj)./T2vals(nphantom))) .* exp(-b_value .* (diffusion_gradient(ii,:) * D(:,:,nphantom) * diffusion_gradient(ii,:)')).*phase(nphantom);
            else
                I(:,:,ii,jj)=I(:,:,ii,jj)+(P{nphantom} .* exp(-T2prep(jj)./T2vals(nphantom))) .* exp(-b_value .* (diffusion_gradient(ii,:) * D(:,:,nphantom) * diffusion_gradient(ii,:)'));
                
            end
        end
        
    end
end

if visualize == 1
% TO DO ...
end

end

