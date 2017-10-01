%
function I=VFA_TSE_T2_T1_phantom(N,T2vals,T1vals,T2prep, TEes, ETL, varargin )
%N : resolution
% T2vals = vector of T2 values
% T1vals = vecotr of T1 values

% TE: vector of TE
% echo_time: vector of echo time values 

%optional: visualize true/false
%% EXAMPLE 
%I=VFA_TSE_T2_T1_phantom(128,[20 40 60 80 100].*1e-3,[300 500 800 1000 1500].*1e-3,[10:3:97].*0.001,4.*0.001,30,1);

%%
assert(length(T2vals)==length(T1vals)); %want same number of bvals as T1vals (= number of phantoms)

if nargin>
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
    
    
    TSE_Mxy_modulation(nphantom,:) = CUBE_Mxy_calculation_lite(T1vals(nphantom)*1e3, T2vals(nphantom)*1e3, TEes*1e3, ETL); % this function take time unit of ms instead of s
end




% make simulation images
I=zeros(N,N,ETL,length(T2prep));
if complexsim
  for nphantom=1:length(T2vals)
    phase(nphantom)=exp(1i*rand*2*pi);
    fprintf('phase phantom %d is %4.2f \n',nphantom,angle(phase(nphantom)))
  end  
end

for ii=1:ETL
    for jj=1:length(T2prep)
        
        for nphantom=1:length(T2vals) % for all separate phantoms
            
            if complexsim %add random phase to phantom
                I(:,:,ii,jj)=I(:,:,ii,jj)+(P{nphantom} .* exp(-T2prep(jj)./T2vals(nphantom))) .* TSE_Mxy_modulation(nphantom,ii).*phase(nphantom);
            else
                I(:,:,ii,jj)=I(:,:,ii,jj)+(P{nphantom} .* exp(-T2prep(jj)./T2vals(nphantom))) .* TSE_Mxy_modulation(nphantom,ii);
            end
        end
        
    end
end

if visualize == 1
% TO DO ...
end

end

