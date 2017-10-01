function D= sortArray(MR)
% function to sort MR.Data if normal MR.SortData crashes
% sorts based on dynamics, averages and channels... (other things would not be sorted
% correctly as of now) 


ky=MR.Parameter.Labels.Index.ky((MR.Parameter.Labels.Index.typ==1));
kz=MR.Parameter.Labels.Index.kz(MR.Parameter.Labels.Index.typ==1);
chan=MR.Parameter.Labels.Index.chan(MR.Parameter.Labels.Index.typ==1);
param1=MR.Parameter.Labels.Index.dyn(MR.Parameter.Labels.Index.typ==1);
param2=MR.Parameter.Labels.Index.aver(MR.Parameter.Labels.Index.typ==1);

MRD=MR.Data;

assert(numel(chan)==length(MRD),'subsets of data not yet supported ')


nkx=size(MR.Data,1);

minky=min(ky);
maxky=max(ky);
nky=length(minky:maxky);
kyunique=unique(ky);

minkz=min(kz);
maxkz=max(kz);
nkz=length(minkz:maxkz);
kzunique=unique(kz);

minchan=min(chan);
maxchan=max(chan);
chanunique=unique(chan);
nchans=length(chanunique);

p1unique=unique(param1+1);
minp1=min(p1unique);
maxp1=max(p1unique);
np1=length(p1unique);

p2unique=unique(param2+1);
minp2=min(p2unique);
maxp2=max(p2unique);
np2=length(p2unique);

assert(minp2==1) 
assert(minp1==1)

%input chan number, outputs index of chan in inuqie chan list
chantransform=zeros(numel(chanunique),1);
for ii=1:numel(chanunique)
chantransform(chanunique(ii))=ii; %[0,1,2,3,4,]
end

%input chan number, outputs index of chan in inuqie chan list
p1indextransform=zeros(numel(p1unique),1);
for ii=1:numel(p1unique)
p1indextransform(p1unique(ii))=ii; %[0,1,2,3,4,]
end


%input chan number, outputs index of chan in inuqie chan list
p2indextransform=zeros(numel(p2unique),1);
for ii=1:numel(p2unique)
p2indextransform(p2unique(ii))=ii; %[0,1,2,3,4,]
end

D=zeros(nkx,nky,nkz,nchans,np1,np2);

kyindex=ky-minky+1;
kzindex=kz-minkz+1;
chanindex= chantransform(chan); % should be able to deal with eg chan=2,3,5,6
np1index= p1indextransform(param1+1); % TO DO: CHANGE: 
np2index=p2indextransform(param2+1); 

fprintf('sorting MR.Data \n')
fprintf('number of unique channels %d \n',nchans)
fprintf('number of unique dynamics %d \n',np1)
fprintf('number of unique averages %d \n',np2)

tic
for ii=1:size(MRD,2)
    if mod(ii,1e3)==0;waitbar(ii/size(MRD,2));end
    D(:,kyindex(ii),kzindex(ii),chanindex(ii),np1index(ii),np2index(ii))=MR.Data(:,ii);
end
toc

end