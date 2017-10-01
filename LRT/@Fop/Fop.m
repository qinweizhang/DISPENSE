function  res = Fop(imsize)
% only for 2D
res.adjoint = 0;
res.imsize=imsize; 
res = class(res,'Fop');

