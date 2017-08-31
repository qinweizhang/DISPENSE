function rand_phase = random_phase_map(kx_dim, ky_dim, kz_dim,width)

%width bigger = smoother

N = [kx_dim, ky_dim, kz_dim]; % size in pixels of image
F = 3;        % frequency-filter width
[X,Y, Z] = ndgrid(1:N(1),1:N(2),1:N(3));
p = min(X-1,N(1)-X+1);
q = min(Y-1,N(2)-Y+1);
r = min(Z-1,N(3)-Z+1);
H = exp(-width*(p.^2+q.^2+r.^2)/F^2);
rand_phase = real(ifftn(H.*fftn(randn(N))));

end