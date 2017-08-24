function b = is2Darray(x)

if length(size(x)) == 2
    if (size(x,1) ~= 1 && size(x,2) ~= 1)
        b = 1;
    else
        b = 0;
    end
else
    b = 0;
end