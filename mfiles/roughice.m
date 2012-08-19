% ROUGHICE   Demonstrate that the SIA code SIAFLAT can deal adaptively
% with terrible surface topography.

% space-time grid parameters
J = 60;  K = J;
Lx = 500e3;  Ly = Lx;  dx = 2 * Lx / J;  dy = 2 * Ly / K;
secpera = 31556926;
x = linspace(-Lx,Lx,J+1);  y = x;
[xx,yy] = meshgrid(x,y);

% construct and plot the worst possible ice sheet
H0 = 3000 + 1000 * rand(J+1,K+1);
H0(xx.^2 + yy.^2 > 300000^2) = 0;
figure(1),  surf(x/1000,y/1000,H0)
xlabel('x  (km)'), ylabel('y  (km)')%, zlabel('h(x,y)   (m)')
title('initial surface')
%print -dpng roughinitial.png

% run SIA code
fprintf('initial ice volume:   %.6e km^3\n',sum(sum(H0))*dx*dy/1e9)
tfyears = 50;
[H,dtlist] = siaflat(Lx,Ly,J,K,H0,0.2*secpera,tfyears*secpera);
fprintf('final ice volume (at t = %.2f a):   %.6e km^3\n',tfyears,sum(sum(H))*dx*dy/1e9)

% show final state and adaptive time-stepping
figure(2),  surf(x/1000,y/1000,H)
xlabel('x  (km)'), ylabel('y  (km)')%, zlabel('h(x,y)   (m)')
title(['final surface at final time = ' num2str(tfyears) ' (a)'])
%print -dpng roughfinal.png

figure(3),  plot(dtlist/secpera,'o')
xlabel('n  (time step count)'), ylabel('dt  (a)')
title('adaptive time steps taken')
%print -dpng roughtimesteps.png

