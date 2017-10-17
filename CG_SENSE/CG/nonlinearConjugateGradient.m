function [z, z_iter] = nonlinearConjugateGradient(f,df,z0,tol,maxit,cg_method,verbose_flag,reference_norm)

% usage:
%   [z, z_iter] = nonlinearConjugateGradient(f,df,z0,tol,maxit,cg_method,verbose_flag,reference_norm)
%
% Finds the minimum of the function specified by the function handle "f"
% using a nonlinear conjugate gradient algorithm. This algorithm is based on
% the freely available code from Michael Overton "NLCG 1.0".
%
% f = function handle of the cost function
% df = function handle of the derivative of the cost function
% z0 = set a starting point
% tol = stopping tolerance
% maxit = maximum iterations
% verbose = use 0 to switch to silent mode
% cg_method = set the CG-Method


persistent h;

if nargin<=3 || isempty(tol)
    tol = 1e-6;
end
if nargin<=5 || isempty(cg_method)
    cg_method = 'fr+pr';
end
if nargin<=6 || isempty(verbose_flag)
    verbose_flag = 2;
end
if nargin<=7 || isempty(reference_norm)
    reference_norm = 1;
end


c1 = 0.05;
backtrack_step = 0.8;

% strongwolfe = 1;
% wolfe1 = 0;
% wolfe2 = 0.5; 
frec = nan; % so defined in case of immediate return
alpharec = nan;
z = z0;

fz = f(z);
dfz = df(z);

if size(dfz,2) > size(dfz,1) % error return is appropriate here
    error('gradient must be returned as a column vector, not a row vector');
end
gnorm = l2norm(dfz);
if fz == inf % better not to generate an error return
    if verbose_flag > 0
        fprintf('nlcg: f is infinite at initial iterate\n');
    end
    return
elseif isnan(fz)
    if verbose_flag > 0
        fprintf('nlcg: f is nan at initial iterate\n');
    end
    return
elseif gnorm < tol
    if verbose_flag > 0
       fprintf('nlcg: tolerance on gradient satisfied at initial iterate\n');
    end
    return
end

p = -dfz;  % start with steepest descent

if nargout==2
    z_iter = cell(1,maxit);
end
for iter = 1:maxit
    
    gtp = real(dotprod(dfz,p));
    if  gtp >= 0 || isnan(gtp)
        if verbose_flag > 0
            fprintf('Not descent direction, taking gradient descent at iteration %d, f(z) = %1.5e, gnorm = %1.5e\n', iter, fz, gnorm);
        end
        p = -dfz;
    end
    
    gprev = dfz;

    % backtracking line search
    alpha = 1;
    if (fz + c1*alpha*gtp) > 0
        fy = f(z,alpha,p,1);
        counter = 1;
        while ( (fy > (fz + c1*alpha*gtp)) ) % sufficient decrease condition aka first Wolfe condition or Armijo rule, second Wolfe condition is too expensive
            alpha = backtrack_step * alpha;
            fy_old = fy;
            fy = f(z,alpha,p,0);
            counter = counter + 1;
            if counter>=200
                break;
                alpha = 0;
            end
        end
        
    else % rare case
        fy_old = Inf;
        fy = f(z,alpha,p,1);
        counter = 1;
        while fy < fy_old
            alpha = backtrack_step * alpha;
            fy_old = fy;
            fy = f(z,alpha,p,0);
            counter = counter + 1;
            if counter>=200
                break;
                alpha = 0;
            end
        end
        alpha = alpha / backtrack_step;
        
    end
    
    z = z + alpha*p;
    
    fz = fy;
    dfz = df(z);
    
    
% original line search method by Michael Overton (backtracking seams to work as well and is faster)
%     if strongwolfe % strong Wolfe line search is usually essential (default)
%         [alpha, z, fz, g, fail] = ...
%             linesch_sw(z, fz, g, p, pars, wolfe1, wolfe2, fvalquit, verbose_flag);
%     else  % still, leave weak Wolfe as an option for comparison
%         [alpha, z, fz, g, fail] = ...
%             linesch_ww(z, fz, g, p, pars, wolfe1, wolfe2, fvalquit, verbose_flag);
%     end


    gnorm = l2norm(dfz);
    frec(iter) = fz;
    alpharec(iter) = alpha;
    
    if verbose_flag==1
        if isempty(h) || ~ishandle(h) || ~strcmp(get(h,'Tag'),'nlcg_figure');
            h = figure('Tag','nlcg_figure');
        end
        if gcf~=h
            set(0,'CurrentFigure',h);
        end
        if isvector(z)
            plot(abs(z));
        elseif is2Darray(z)
            imagesc(abs(z));
            colormap gray;
        elseif length(size(z))==3
            imagesc(abs(array2mosaic(z)));
            colormap gray;
        end
        axis off;
        title(sprintf('iter %d: step = %1.5e, f(z) = %1.5e, gnorm = %1.5e\n', iter, alpha*gnorm/reference_norm, fz, gnorm));
        drawnow;
        fprintf('iter %d: step = %1.5e, f(z) = %1.5e, gnorm = %1.5e\n', iter, alpha*gnorm/reference_norm, fz, gnorm);
        
        figure(3000);
        if isvector(z)
            plot(angle(z));
        elseif is2Darray(z)
            imagesc(angle(z),[-pi pi]); colormap jet
            colormap gray;
        elseif length(size(z))==3
            imagesc(angle(array2mosaic(z)),[-pi pi]); colormap jet
            colormap gray;
        end
        axis off;
        title(sprintf('iter %d: step = %1.5e, f(z) = %1.5e, gnorm = %1.5e\n', iter, alpha*gnorm/reference_norm, fz, gnorm));
        drawnow;
        fprintf('iter %d: step = %1.5e, f(z) = %1.5e, gnorm = %1.5e\n', iter, alpha*gnorm/reference_norm, fz, gnorm);

        
    elseif verbose_flag==2
        fprintf('iter %d: step = %1.5e, f(z) = %1.5e, gnorm = %1.5e\n', iter, alpha*gnorm/reference_norm, fz, gnorm);
    end
    
    if (alpha*gnorm/reference_norm) <= tol
        if verbose_flag > 0
            fprintf('step length below tolerance, quit at iteration %d, f(z) = %1.5e\n', iter, fz);
        end
        return
    end
    
    y = dfz - gprev;
    
    if strcmp(cg_method,'pr')
        nmgprevsq = dotprod(gprev,gprev);
        beta = dotprod(dfz,y)/nmgprevsq;  % Polak-Ribiere-Polyak
    elseif strcmp(cg_method,'fr')
        nmgprevsq = dotprod(gprev,gprev);
        beta = dotprod(dfz,dfz)/nmgprevsq;  % Fletcher-Reeves  
    elseif strcmp(cg_method,'fr+pr') % combined PR-FR (suggested by Gilbert-Nocedal)
    % ensures beta <= |beta_fr|, allowing proof of % global convergence,
    % but avoids inefficiency of FR which happens when beta_fr gets stuck near 1
        nmgprevsq = dotprod(gprev,gprev);
        beta_pr = dotprod(dfz,y)/nmgprevsq;  % Polak-Ribiere-Polyak
        beta_fr = dotprod(dfz,dfz)/nmgprevsq;  % Fletcher-Reeves  
        if beta_pr < -beta_fr  
            if verbose_flag > 1
                fprintf('*** truncating beta_pr = %1.5e to -beta_fr = %1.5e\n', beta_pr, -beta_fr);
            end
            beta = -beta_fr;  
        elseif beta_pr > beta_fr
            if verbose_flag > 1
                fprintf('*** truncating beta_pr = %1.5e to +beta_fr = %1.5e\n', beta_pr, beta_fr);
            end
            beta = beta_fr;
        else
            beta = beta_pr;
        end
    elseif strcmp(cg_method,'hs') % Hestenes-Stiefel
        beta = dotprod(dfz,y)/dotprod(p,y);
    elseif strcmp(cg_method,'dy') % Dai-Yuan
        beta = dotprod(dfz,dfz)/dotprod(p,y);
    elseif strcmp(cg_method,'hz') % Hager-Zhang
        pty = dotprod(p,y); % p is called d in their paper
        theta = 2*dotprod(y,y)/pty;
        beta_hz = dotprod(y-theta*p,dfz)/pty;
        eta = -1/(l2norm(p)*min(0.01,l2norm(gprev)));
        beta = max(beta_hz, eta);
    elseif strcmp(cg_method,'-') % Steepest Descent
        beta = 0;
    else
        error('No such cg_method');
    end
    
    p = beta*p - dfz;

    if nargout==2
        z_iter{iter} = z;
    end
    
end % for loop
if verbose_flag > 0
    fprintf('%d iterations reached, f(z) = %1.5e, gnorm = %1.5e\n', maxit, fz, gnorm);
end