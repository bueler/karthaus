function T = heatadapt(D,J,K,tf)
% HEATADAPT  Adaptive time-stepping explicit method for heat equation
%   T_t = D (T_xx + T_yy)
% on -1 < x < 1, -1 < y < 1, for 0 < t < tf
% Compare HEAT.  This version uses the same gaussian initial condition.
% Calling sequence slightly different because time step dt is computed
% and is not an input.
% Usage:
%   T = heatadapt(D,J,K,tf)
% where
%   T   = approximate solution at tf
%   D   = diffusivity coeff
%   J,K = number of points in x,y directions, resp.
%   tf  = final time
% Examples:  Compare these
%   >> heatadapt(1.0,50,50,0.05);
%   >> heat(1.0,50,50,0.00025,200);

% spatial grid and initial condition:
dx = 2 / J;    dy = 2 / K;
[x,y] = ndgrid(-1:dx:1, -1:dy:1); % (J+1) x (K+1) grid in x,y plane
T = exp(-30*(x.*x + y.*y));

fprintf('  doing explicit steps adaptively on 0.0 < t < %.3f\n',tf)
t = 0.0;    count = 0;
while t < tf
   % stability requires  1 - 2 nu_x - 2 nu_y >= 0,  which is the
   %   following as a time-step restriction
   dt0 = 0.25 * min(dx,dy)^2 / D;
   dt = min(dt0, tf - t);  % do not go past tf
   nu_x = dt * D / (dx*dx);    nu_y = dt * D / (dy*dy);
   T(2:J,2:K) = T(2:J,2:K) + ...
       nu_x * ( T(3:J+1,2:K) - 2 * T(2:J,2:K) + T(1:J-1,2:K) ) + ...
       nu_y * ( T(2:J,3:K+1) - 2 * T(2:J,2:K) + T(2:J,1:K-1) );
   t = t + dt;
   count = count + 1;
   fprintf('.')
end
fprintf('\n  completed N = %d steps, average dt = %.7f\n',count,tf/count)

surf(x,y,T),  shading('interp'),  xlabel x,  ylabel y
