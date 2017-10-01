%
function I=diffusion_T1_phantom(N,ADC,T1vals,TE,b,varargin)
%N : resolution
% ADC = vector of bvalues
% T1vals = vecotr of T1 values

% TE: vector of TE
% b: vector of b values 

%optional: visualize true/false
%% EXAMPLE 
% I=diffusion_T1_phantom(100,[0.1 0.5 1],[200, 400 600],[1,2,3,4,5].*0.001,[100,200,300],1)

%%
assert(length(ADC)==length(T1vals)); %want same number of bvals as T1vals (= number of phantoms)

if nargin>5
    visualize = varargin{1};
else
    visualize =false;
end

%% START CODE
for nphantom=1:length(ADC)
    E(1)=0.4+0.6*rand; %intensity
    E(2)= 0.1+rand*0.2; %length
    E(3)= 0.1+rand*0.2; %width
    E(4)=-0.5+(rand); %x-coord of middle 
    E(5)=-0.5+(rand); %y-coord of middle
    E(6)= rand*360; %angle   
    P{nphantom} = phantom(E,N) ;
end

% make simulation images
I=zeros(N,N,length(b),length(TE));
for ii=1:length(b)
    for jj=1:length(TE)
        
        for nphantom=1:length(ADC) % for all separate phantoms
            I(:,:,ii,jj)=I(:,:,ii,jj)+(P{nphantom}.*exp(-TE(jj)./T1vals(nphantom))).*exp(-ADC(nphantom).*b(ii));
        end
        
    end
end

if visualize == 1
% TO DO ...
end

end

