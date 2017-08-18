function  norms = colnorm2(x)
%Calculate l2-norm for every columns in a 2D matrix


    assert(length(size(x))<=2);
    
    for col = 1:size(x,2)
        norms(col)  = norm(x(:,col),2);
    end

end