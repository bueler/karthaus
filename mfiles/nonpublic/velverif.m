function velverif(J)

FIXME: TO DO

% grid
L = 800e3;  dx = L/J;  dy = dx; 
[x,y] = meshgrid(-L:dx:L,-L:dy:L);

% constants as in halfar.m:
H0 = 3600;  R0 = 750e3;  secpera = 31556926;
n = 3;  alpha = 1/9;  beta = 1/18;
Gamma = 9.0177e-13;  t0 = (beta/Gamma) * (7/4)^3 * (R0^4/H0^7);

% 
r = sqrt(x.*x + y.*y);
r = r / R0;  t = t / t0;
inside = max(0, 1 - (r / t^beta).^((n+1) / n));
H = H0 * inside.^(n / (2*n+1)) / t^alpha;


csurf = csurf * secpera;  % units of m/a
csurf(csurf < 1e-4) = 1e-4;  % cut off below at 0.1 mm/a so log10() works better
imagesc(x(2:J)/1000,y(2:J)/1000,log10(flipud(csurf))), axis square, colorbar
xlabel x, ylabel y, title('surface velocity (log10 of m/a) for initial geometry')
%print -dpng antinitcsurf.png

