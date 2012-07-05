function ant(doplot,E)
% ANT  simulate Antarctic ice sheet flow using buildant.m
%   to extract data from re-gridded SeaRISE-Antarctic data
%   (see Ant50km.nc for metadata), and siageneral.m to solve SIA
% calls:  buildant.m, siageneral.m, diffusion.m

if nargin < 1, doplot = 1; end
if nargin < 2, E = 3; end  % default enhancement factor

[x,y,lat,lon,prcp,thk,topg,usrf] = buildant(0);  % read input data from NetCDF; no plot

% grid info
Lx = (max(x) - min(x)) / 2;    Ly = (max(y) - min(y)) / 2;
J = length(x) - 1;    K = length(y) - 1;
dx = 2 * Lx / J;    dy = 2 * Ly / K;

fprintf('summary of input data:\n')
fprintf('  dx = %.3f km,  dy = %.3f km\n',dx/1000.0,dy/1000.0)
fprintf('  thickness     [min,max] = [%8.2f,%8.2f] m\n',min(min(thk)),max(max(thk)))
fprintf('  bed elevation [min,max] = [%8.2f,%8.2f] m\n',min(min(topg)),max(max(topg)))
fprintf('  precipitation [min,max] = [%8.5f,%8.5f] m a-1\n',min(min(prcp)),max(max(prcp)))

% run-time and time-step (in years)
deltat = 1.0;
tf = 500.0;
NN = 80;  % number of blocks of length tf
fprintf('doing run of %.3fka total, in blocks of %.3fa,\n',...
        NN*tf/1000.0,tf)
fprintf('  with max time-steps of %.3fa ...\n  running:\n',deltat)

% fix units and parameters for actual run
secpera = 31556926;    rho = 910.0;    rhow = 1028.0;
deltat = deltat * secpera;  tf = tf * secpera;  % convert to seconds
M = prcp / secpera;  % alternatively: M = zeros(size(thk));
if any(any(thk<0)), error('negative input thicknesses detected'), end
A = E * 1.0e-16 / secpera;

% solve SIA by doing blocks of tf time, and reporting volume
H = thk;
vol = printvolume(0.0,dx,dy,H);
hinit = getsurface(H,topg,rho,rhow);
for k = 1:NN
  [H,hfinal,dtlist] = siageneral(Lx,Ly,J,K,H,deltat,tf,topg,M,A);
  if any(any(H<0)), error('negative output thicknesses detected'), end
  vol = [vol printvolume(k*tf/secpera,dx,dy,H)];
end

if doplot==0, return; end

figure(1)
imagesc(x/1000,y/1000,flipud(hinit),[0, 4000]), colorbar, axis square
xlabel('x  (km)','fontsize',14), ylabel('y  (km)','fontsize',14)
title('initial surface elevation')
%print -dpdf antinitial.pdf

figure(2)
imagesc(x/1000,y/1000,flipud(hfinal),[0, 4000]), axis square, colorbar
xlabel('x  (km)','fontsize',14), ylabel('y  (km)','fontsize',14)
title('final surface elevation')
%print -dpdf antfinal.pdf

figure(3)
imagesc(x/1000,y/1000,flipud(H-thk),[-1000, 1000]), axis square, colorbar
xlabel x, ylabel y, title('thickness change')

figure(4)
plot((0:NN)*tf/secpera,vol/(1.0e6*1.0e9),'o-','markersize',11,'linewidth',2)
xlabel('t  (a)','fontsize',14), ylabel('volume  (10^6 km^3)','fontsize',14)
grid on
%print -dpdf antvol.pdf

  function vol = printvolume(time,dx,dy,thk)
    vol = sum(sum(thk)) * dx * dy;
    fprintf('  ice volume at time %7.3fka = %.4e km^3\n',time/1000.0,vol/1.0e9)

