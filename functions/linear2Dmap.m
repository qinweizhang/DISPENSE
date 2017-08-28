%Generate 2D linear map
% 
% INPUT
% 
% theta:  gradient direction in 2D plane [-pi pi]
% grad:   gradient strength in direction theta
% xdim:   x dimention of output
% ydim:   y dimention of output
% 
% OUTPUT
% 
% res:    generated 2D map
% 
% Method: project any vector x to a = (cos(theta), sin(theta)), by xp = inv(a'a)*a'"x;
% then give value based on abs(xp)*grad. This can be easily extended to nD
% 
% 
% (c) q.zhang @AMC Amsterdam 2017



function res = linear2Dmap(theta, grad, xdim, ydim)

res = zeros(xdim, ydim);
a = [cos(theta); sin(theta)];
ata = (a'*a);
for x=1:xdim
    for y = 1:ydim
        res(x, y) = grad * abs(a'*[x;y]./ata);
    end
end


end