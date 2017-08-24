function Q = mtimes(A,B)

if A.has_adjoint
    if A.adjoint==0
        Q = A.fhandle(B);
    else
        Q = A.fhandle_adjoint(B);
    end
else
    Q = A.fhandle(B);
end
