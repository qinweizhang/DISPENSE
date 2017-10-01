function [ksp_scaled,scaling]= scaleksp(ksp,imrecon)

tmp = dimnorm(imrecon, 3); % norm in coil dim?
tmpnorm2 = sort(tmp(:), 'ascend');
p100 = tmpnorm2(end);
p90 = tmpnorm2(round(.9 * length(tmpnorm2)));
p50 = tmpnorm2(round(.5 * length(tmpnorm2)));
if (p100 - p90) < 2 * (p90 - p50)
    scaling = p90;
else
    scaling = p100;
end
fprintf('\nScaling: %f\n\n', scaling);

ksp_scaled = ksp ./ scaling;

end

function A2=dimnorm( A, dim, arg )
%dimnorm Compute norm along direction dim
%   arg: second argument passed to norm function (default: 2)

if nargin < 3
    arg = 2;
end


dims = size(A);
dims2 = dims; dims2(dim) = 1;

ldims = 1:length(dims);

ldims2 = [dim, ldims]; ldims2(dim+1) = [];

B = reshape(permute(A, ldims2), dims(dim), []);
y = vecnorm(B, arg);
A2 = reshape(y, dims2);
end


function [ y ] = vecnorm( A, arg )
%vecnorm Compute norm of each column vector in a matrix
%   arg: second argument passed to norm function (default: 2)

if nargin < 2
    arg = 2;
end

% avoid for loop if possible
if isnumeric(arg) && (abs(arg) ~= Inf)
    y = sum(abs(A).^arg,1).^(1/arg);
else
    [~, ny] = size(A);
    y = zeros(1, ny);
    for ii=1:ny
        a = A(:,ii);
        y(ii) = norm(a,arg);
    end
end

end
