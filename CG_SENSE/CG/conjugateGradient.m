function [z,z_iter] = conjugateGradient(A,b,tol,maxit,z0,verbose_flag)

% function z = conjugateGradient(A,b,tol,maxit,z0)
%
% Solves: A*z=b
% A must be hermitian and positive definite.
%
% tol = tolerance of the method (default: 1e-6)
% maxit = maximum number of iterations (default: 20)
% z0 = initial guess (default: 0)

persistent h;

start_at_zero = 1;

if nargin<=2 || isempty(tol)
    tol = 1e-6;
end
if nargin<=3 || isempty(maxit)
    maxit = 20;
end
if nargin<=4 || isempty(z0)
    if isnumeric(b)
        z = zeros(size(b));
    else
        eval(['z = ',class(b),'(size(b));']);
    end
elseif isnumeric(z0) && isscalar(z0) && z0==0
    z = zeros(size(b));
else
    z = z0;
    start_at_zero = 0;
end
if nargin<=5
    verbose_flag = 1;
end

%nb = sqrt(sum(col(conj(b).*b)));
%tolb = tol * nb;
if start_at_zero
    d = b;
else
    d = b - A*z;
end
r = d;

normrr0 = sum(col(conj(r) .* r));
normrr = normrr0;

for n=1:maxit
    Ad = A*d;
    tmp = conj(d) .* Ad;
    alpha = normrr/sum(col(tmp));
    z = z + alpha*d;
    r = r - alpha*Ad;
    normrr2 = sum(col(conj(r) .* r));
	beta = normrr2/normrr;
    normrr = normrr2;
	d = r + beta*d;
    
    if nargout==2
        z_iter{n} = z;
    end
    
    if verbose_flag==1
        if isempty(h) || ~ishandle(h) || ~strcmp(get(h,'Tag'),'cg_figure');
            h = figure('Tag','cg_figure');
        end
        if gcf~=h
            set(0,'CurrentFigure',h);
        end
        if isnumeric(z) && isvector(z)
            plot(abs(z));
        elseif isnumeric(z) && is2Darray(z)
            imagesc(abs(z));
            colormap gray;
        elseif isnumeric(z) && length(size(z))==3
            imagesc(abs(array2mosaic(z)));
            colormap gray;
        else
            imagesc(z);
            colormap gray;
        end
        axis off;
        title(['Iteration #', num2str(n),': current tolerance = ', num2str(sqrt(normrr/normrr0)), '.']);
        drawnow;
        
         figure(2000);
        if isnumeric(z) && isvector(z)
            plot(angle(z));
        elseif isnumeric(z) && is2Darray(z)
            imagesc(angle(z),[-pi pi]);
            colormap jet;
        elseif isnumeric(z) && length(size(z))==3
            imagesc(angle(array2mosaic(z)),[-pi pi]);
            colormap jet;
        else
            imagesc(z,[-pi pi]);
            colormap jet;
        end
        axis off;
        title(['Iteration #', num2str(n),': current tolerance = ', num2str(sqrt(normrr/normrr0)), '.']);
        drawnow;
        
    elseif verbose_flag==2
        fprintf(['Iteration #', num2str(n),': current tolerance = ', num2str(sqrt(normrr/normrr0)), '\n']);
    end
    
    if sqrt(normrr/normrr0) < tol
        break;
    end
    
end
