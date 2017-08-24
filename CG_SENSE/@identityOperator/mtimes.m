function Q = mtimes(A,B)

if strcmp(class(A),'identityOperator')
    Q = B;
    
% now B is the operator and A is the vector
elseif strcmp(class(B),'identityOperator')
    Q = mtimes(B',A')';
    
end