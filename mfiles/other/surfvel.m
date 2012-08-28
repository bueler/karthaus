function [Urtx,Urty,Uupx,Uupy] = surfvel(x,y,h,thk,n,A)

FIXME:  TO TEST

% SURFVEL  Compute surface velocity from non-sliding SIA with isothermal Glen law:
%           2           n   n+1         n-1
%   U = -  --- A (rho g)   H    |grad h|    grad h
%          n+1
% See formulas (5.84) in Greve & Blatter (2009), in case (v_bx,v_by)=(0,0) and
% when A(T)=A is constant.  Computes x,y components of velocity on the staggered
% grid using Mahaffy scheme.
%
% form:
%   [Urtx,Urty,Uupx,Uupy] = surfvel(x,y,h,thk,n,A)
% where
%   x = J+1 length array of equally-spaced x values
%   y = K+1 length array of equally-spaced y values
%   h = surface elevation = (J+1) x (K+1) array
%   thk = H = ice thickness
%   n = Glen exponent
%   A = ice softness
% all inputs in SI units: x,y,h,thk (m);  A (Pa^-n s^-1)

% note:  Array  h = h(1:J+1,1:K+1)   is size  J+1 x K+1, and includes boundary
%        of rectangle.  Thus interior points are  h(2:J,2:K)  and is size
%        J-1 x K-1.  Similarly the staggered gradient components dhd?_?? are
%        also size  J-1 x K-1.

% grid spacing
if norm(diff(diff(x)))/norm(diff(x)) > 1.0e-8
  error('input array x is not equally-spaced'), end
if norm(diff(diff(y)))/norm(diff(y)) > 1.0e-8
  error('input array y is not equally-spaced'), end
dx = x(2) - x(1);   dy = y(2) - y(1);

% compute surface gradient at staggered grid locations
dhdx_rt = (h(3:J+1,2:K) - h(2:J,2:K)) / dx;
dhdy_rt = (h(2:J,3:K+1) + h(3:J+1,3:K+1) - h(2:J,1:K-1) - h(3:J+1,1:K-1)) / (4*dy);
dhdx_up = (h(3:J+1,2:K) + h(3:J+1,3:K+1) - h(1:J-1,2:K) - h(1:J-1,3:K+1)) / (4*dx);
dhdy_up = (h(2:J,3:K+1) - h(2:J,2:K)) / dy;
gradh_rt = sqrt(dhdx_rt.^2 + dhdy_rt.^2);  % magnitudes
gradh_up = sqrt(dhdx_up.^2 + dhdy_up.^2);

% staggered thicknesses are averages from regular grid
H_rt = 0.5 * (thk(2:J,2:K) + thk(3:J+1,2:K));
H_up = 0.5 * (thk(2:J,2:K) + thk(2:J,3:K+1));

% CC = overall constant in velocity
g = 9.81;   rho = 910.0;
CC = (2 * A * (rho * g)^n) / (n+1);

% result
scalar_rt = CC * gradh_rt.^n .* H_rt.^(n+1);
scalar_up = CC * gradh_up.^n .* H_up.^(n+1);
Urtx = scalar_rt .* dhdx_rt;
Urty = scalar_rt .* dhdy_rt;
Uupx = scalar_up .* dhdx_up;
Uupy = scalar_up .* dhdy_up;


