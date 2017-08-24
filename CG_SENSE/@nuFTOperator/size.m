function [s1,s2] = size(A,n)

t1 = [A.trajectory_length*A.numCoils 1];
t2 = A.imageDim;

if A.adjoint
    tmp = t1;
    t1 = t2;
    t2 = tmp;    
end

if nargin==1
    s1 = t1;
    s2 = t2;
elseif nargin==2
    s2 = [];
    if n==1
        s1 = t1;
    elseif n==2
        s1 = t2;
    end
end
