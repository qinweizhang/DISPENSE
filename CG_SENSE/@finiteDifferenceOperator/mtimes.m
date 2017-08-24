function Q = mtimes(A,B)

if strcmp(class(A),'finiteDifferenceOperator')    
    if A.adjoint==1
        if isvector(B)
            Q = zeros(size(B));
            Q(end+1) = 0;
            Q(1) = -B(1);
            Q(2:end-1) = -diff(B);
            Q(end) = B(end);
            
        else
            s = [1:length(size(B))];
            s(A.direction) = [];
            s = [A.direction, s];
            B = permute(B, s);

            ddots = '';
            for k=1:length(size(B))-1
                ddots = [',:', ddots];
            end

            E1 = eval(['B(1', ddots, ');']);
            Eend = eval(['B(end', ddots, ');']);

            Q = cat(1, -E1, -diff(B,1,1), Eend);
            Q = ipermute(Q, s);
        end

    else
        if isvector(B)
            Q = diff(B);
        else
            Q = diff(B,1,A.direction);
        end
        
    end

% now B is the operator and A is the vector
else
    Q = mtimes(B',A')';
    
end