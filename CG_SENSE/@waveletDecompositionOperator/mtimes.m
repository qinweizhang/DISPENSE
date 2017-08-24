function Q = mtimes(A,B)

if strcmp(class(A),'waveletDecompositionOperator')
    
    if A.adjoint==1
        if A.space==1
            Q = waverecX(B,A.params.sizes,A.wname);
            
        elseif A.space==2
            Q = waverec2X(B,A.params.sizes,A.wname);
            
        elseif A.space==3
            X.sizeINI = A.params.sizeINI;
            X.level = A.params.level;
            X.mode = A.params.mode;
            X.filters = A.params.filters;
            X.sizes = A.params.sizes;
            
            ps = prod(A.sizes_all,2);
            Y = cell(1,length(ps));
            startpt = 1;
            for k=1:length(ps)
                Y{k} = reshape(B(startpt:startpt+ps(k)-1),A.sizes_all(k,:));
                startpt = startpt + ps(k);
            end
            X.dec = Y;
            Q = waverec3X(X);
            
        end
            
        
    else
        
        if A.space==1
            Q = wavedecX(B,A.N,A.wname);
            
        elseif A.space==2
            Q = wavedec2X(B,A.N,A.wname);
            
        elseif A.space==3
            X = wavedec3X(B,A.params.level,A.wname);
            X = X.dec;
            
            ps = prod(A.sizes_all,2);
            Q = zeros(sum(ps),1);
            startpt = 1;
            for k=1:length(X)
                Q(startpt:startpt+ps(k)-1) = col(X{k});
                startpt = startpt + ps(k);
            end
            
        end
                
    end

% now B is the operator and A is the vector
else
    Q = mtimes(B',A')';
    
end