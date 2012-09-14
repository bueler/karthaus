% ASHTRAY  Create an ice cap which sits on an idealized bedrock
% topography which is like a cup, and has a trough in one direction.
% So it looks like an ashtray.  Use SIAGENERAL to evolve this ice cap
% from zero ice to near steady state.  This is something to model both
% as a region and as a whole ice sheet.  This is for idea 2 in ideas.txt.

% space-time grid parameters
J = 100;  K = J;
Lx = 900e3;  Ly = Lx;
dx = 2 * Lx / J;  dy = 2 * Ly / K;
secpera = 31556926;
x = linspace(-Lx,Lx,J+1);
y = linspace(-Ly,Ly,K+1);
[xx,yy] = meshgrid(x,y);

% construct bedrock and initial thickness
L = 750e3;
H0 = zeros(size(xx));
rr = sqrt(xx.^2+yy.^2);
b = - 500 * cos(1.2 * pi * rr / L) + 1000;
trough = max(0,xx) * (600/Lx) .* exp(-(3*yy/L).^2);
b = b - trough;
M = (0.3/secpera) * (1 - 2 * (rr / L).^2);

% show what we have so far
%figure(1),  surf(x/1000,y/1000,b)
%xlabel('x  (km)'), ylabel('y  (km)')
%title('bedrock elevation b(x,y)  (m)')
%figure(2),  surf(x/1000,y/1000,M*secpera)
%xlabel('x  (km)'), ylabel('y  (km)')
%title('mass balance M(x,y)  (m/a)')

% run SIA code
%fprintf('initial ice volume:   %.6e km^3\n',sum(sum(H0))*dx*dy/1e9)
E = 5;
A = E * 1.0e-16 / secpera;
tf = 40000.0 * secpera;
deltat = 10.0 * secpera;
[H,hfinal,dtlist] = siageneral(Lx,Ly,J,K,zeros(size(b)),deltat,tf,b,M,A);
fprintf('minimum time step = %.6f a,  maximum time step = %.6f a\n',...
   min(dtlist)/secpera,max(dtlist)/secpera)
fprintf('final ice volume (at t = %.2f a):   %.6e km^3\n',tf/secpera,sum(sum(H))*dx*dy/1e9)

% show result and time steps
figure(3),  contour(x/1000,y/1000,b+H), view(2), axis equal
xlabel('x  (km)'), ylabel('y  (km)')
figure(4), plot((deltat:deltat:tf)/secpera,dtlist/secpera)
xlabel('t  (a)'), title('time step  (a)')

