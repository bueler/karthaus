function [H,dtlist] = siaflat(Lx,Ly,J,K,H0,deltat,tf)
% SIAFLAT  numerical solution of isothermal n=3 SIA with flat bed and
%          no accumulation; uses Mahaffy (1976) method to evaluate map-plane
%          diffusivity
% equation:
%   H_t = div (D grad H)
% where
%   H = ice thickness = ice surface elevation
%   D = Gamma H^5 |grad H|^2 = solution-dependent (nonlinear) diffusivity
%   Gamma = 2 A (rho g)^3 / 5
% form
%   [H,dtlist] = siaflat(Lx,Ly,J,K,H0,deltat,tf)
% where
%   Lx,Lx = half lengths of rectangle in x,y directions
%   J,K = number of subintervals in x,y directions
%   H0 = initial thickness, a (J+1)x(K+1) array
%   deltat = major time step
%   tf = final time
% outputs
%   H = numerical approx of thickness at final time
%   dtlist = list of time steps used adaptively in diffusion.m
% calls diffusion.m, which does adaptive explicit time-stepping
%   within the major time step
% example: see verifysia.m

% constants
g = 9.81;    rho = 910.0;    secpera = 31556926;
A = 1.0e-16/secpera;    Gamma  = 2 * A * (rho * g)^3 / 5; % see Bueler et al (2005)
H = H0;

dx = 2 * Lx / J;   dy = 2 * Ly / K;
N = ceil(tf / deltat);    deltat = tf / N;
j  = 2:J;    k = 2:K;  % interior indices
nk = 3:K+1;    sk = 1:K-1;    ej = 3:J+1;    wj = 1:J-1; % north,south,east,west

fprintf('SIA solver calling diffusion(), computing diffusivity by Mahaffy ...\n')
t = 0;   dtlist = [];
for n=1:N
  % staggered grid thicknesses
  Hup = 0.5 * ( H(j,nk) + H(j,k) );   Hdn = 0.5 * ( H(j,k) + H(j,sk) ); % up and down
  Hrt = 0.5 * ( H(ej,k) + H(j,k) );   Hlt = 0.5 * ( H(j,k) + H(wj,k) ); % right and left
  % staggered grid value of |grad h|^2 = "alpha^2"
  a2up = (H(ej,nk) + H(ej,k) - H(wj,nk) - H(wj,k)).^2 / (4*dx)^2 + ...
         (H(j,nk) - H(j,k)).^2 / dy^2;
  a2dn = (H(ej,k) + H(ej,sk) - H(wj,k) - H(wj,sk)).^2 / (4*dx)^2 + ...
         (H(j,k) - H(j,sk)).^2 / dy^2;
  a2rt = (H(ej,k) - H(j,k)).^2 / dx^2 + ...
         (H(ej,nk) + H(j,nk) - H(ej,sk) - H(j,sk)).^2 / (4*dy)^2;
  a2lt = (H(j,k) - H(wj,k)).^2 / dx^2 + ...
         (H(wj,nk) + H(j,nk) - H(wj,sk) - H(j,sk)).^2 / (4*dy)^2;     
  % Mahaffy evaluation of staggered grid diffusivity
  %   D = Gamma H^{n+2} |grad h|^{n-1}
  Dup = Gamma * Hup.^5 .* a2up;  Ddn = Gamma * Hdn.^5 .* a2dn;
  Drt = Gamma * Hrt.^5 .* a2rt;  Dlt = Gamma * Hlt.^5 .* a2lt;
  % call *adaptive* diffusion() to time step H
  [H,dtadapt] = diffusion(Lx,Ly,J,K,Dup,Ddn,Drt,Dlt,H,deltat);
  t = t + deltat;    dtlist = [dtlist dtadapt];
  fprintf('.')
end
fprintf('\nSIA solver done.\n')

% to plot the diffusivity for the final state:
%x = linspace(-Lx+dx,Lx-dx,J-1);  y = linspace(-Ly+dy,Ly-dy,K-1);
%figure(99), surf(x,y,0.25*(Dup+Ddn+Drt+Dlt))

