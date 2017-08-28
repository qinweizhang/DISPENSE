function B = ctranspose(A)

B = A;
if B.adjoint==0
    B.adjoint = 1;
else
    B.adjoint = 0;
end