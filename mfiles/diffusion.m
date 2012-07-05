function [U,dtav] = diffusion(Lx,Ly,J,K,Dup,Ddown,Dright,Dleft,U0,tf,b)
% DIFFUSION  adaptive explicit method for non-constant diffusion equation
%   u_t = div (D grad u)
% form
%   U = diffusion(Lx,Ly,J,K,Dup,Ddown,Dright,Dleft,U0,tf)
% where Lx,Ly are half-widths of rectangular domain
%       D*    are (J-1) x (K-1) matrices with diffusivities for "staggered" grid
%       U0    is (J+1) x (K+1) matrix with initial values on regular grid
% note: no error checking on sizes of D* or U0 matrices
% note: input diffusivities could be time-dependent, but are time-independent
%       here in this simplified implementation; call-back to a diffusivity
%       function would be one way to implement time-dependence D(t,x,y)
% example: compare this result to heatadapt():
%   >> J=50; K=50; D=ones(J-1,K-1);
%   >> [x,y]=ndgrid(-1:2/J:1,-1:2/K:1);
%   >> U0 = exp(-30*(x.*x + y.*y));
%   >> U = diffusion(1.0,1.0,J,K,D,D,D,D,U0,0.05);
%   >> surf(x,y,U), shading('interp'), xlabel x, ylabel y, zlabel u
% major examples: see siaflat.m

% spatial grid and initial condition
dx = 2 * Lx / J;    dy = 2 * Ly / K;
[x,y] = ndgrid(-Lx:dx:Lx, -Ly:dy:Ly); % (J+1) x (K+1) grid in x,y plane
U = U0;
if nargin < 11, b = zeros(size(U0)); end  % non-required argument b will allow
                                          % use for nonflat-bed SIA case

t = 0.0;    count = 0;
while t < tf
   % stability condition as a time-step restriction:
   Dregular = max(max(Dup,Ddown),max(Dleft,Dright));  % array on regular grid
   maxD = max(max(Dregular));  % scalar maximum of D
   dt0 = 0.25 * min(dx,dy)^2 / maxD;
   dt = min(dt0, tf - t);  % do not go past tf
   mu_x = dt / (dx*dx);    mu_y = dt / (dy*dy);
   Ushift = U + b;  % differs from U only if b is not zero
   U(2:J,2:K) = U(2:J,2:K) + ...
       mu_y * Dup    .* ( Ushift(2:J,3:K+1) - Ushift(2:J,2:K) ) - ...
       mu_y * Ddown  .* ( Ushift(2:J,2:K) - Ushift(2:J,1:K-1) ) + ...
       mu_x * Dright .* ( Ushift(3:J+1,2:K) - Ushift(2:J,2:K) ) - ...
       mu_x * Dleft  .* ( Ushift(2:J,2:K) - Ushift(1:J-1,2:K) );
   t = t + dt;    count = count + 1;
end
dtav = tf / count;

