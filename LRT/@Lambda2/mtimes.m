function res = mtimes(a,b)
F=a.F;
G=a.G;
Phi=a.Phi;

if a.adjoint %L'(d) ===> C
    res=G'*((F'*(b))*Phi'); % not totally sure about this operator.
else %L*C= F G C Phi
    res = F*(G*b*Phi);
    res=res.*a.mask;
end