function F = ctranspose(F)
val=~F.adjoint; 
set_MCFop_adjoint(F,val)
% % adj = 1;
% set_MCFop_adjoint(a,adj);
% res = a;

