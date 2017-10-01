% TEMP
function R=loopfunc(AA,Phi,mask,C)
% loop over all time points (can be done in parallel) 
parfor i=1:size(Phi,2)
   Q=Phi(:,i)*Phi(:,i).';
   kpoints_time=bsxfun(@times,AA,mask(:,i));
   result=Q*kpoints_time.';
   results=result.';
   results=results*C.';
   R(:,i)=results; 
end
R=sum(R,2); %??
end