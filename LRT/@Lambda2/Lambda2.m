function  res = Lambda2(F,G,Phi,mask)
% only for 2D
res.F=F;
res.G=G;
res.Phi=Phi;
res.adjoint = 0;
res.mask=mask;
res = class(res,'Lambda2');

