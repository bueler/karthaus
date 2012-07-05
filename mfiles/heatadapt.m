function U = heatadapt(D,J,K,tf)
% adaptive time-stepping explicit method for heat equation
%   u_t = D (u_xx + u_yy)
% on -1 < x < 1, -1 < y < 1, for 0 < t < tf
% uses specific gaussian initial condition
% form as in heat.m:
%   U = heatadapt(D,J,K,tf)
% compare to heat.m:  >> heatadapt(1.0,50,50,0.05);
%                     >> heat(1.0,50,50,0.00025,200);

% spatial grid and initial condition:
dx = 2 / J;    dy = 2 / K;
[x,y] = ndgrid(-1:dx:1, -1:dy:1); % (J+1) x (K+1) grid in x,y plane
U = exp(-30*(x.*x + y.*y));

fprintf('  doing explicit steps adaptively ...\n')
t = 0.0;    count = 0;
while t < tf
   % stability requires  1 - 2 nu_x - 2 nu_y >= 0,  which is the following as
   %   a time-step restriction:
   dt0 = 0.25 * min(dx,dy)^2 / D;
   dt = min(dt0, tf - t);  % do not go past tf
   nu_x = dt * D / (dx*dx);    nu_y = dt * D / (dy*dy);
   U(2:J,2:K) = U(2:J,2:K) + ...
       nu_x * ( U(3:J+1,2:K) - 2 * U(2:J,2:K) + U(1:J-1,2:K) ) + ...
       nu_y * ( U(2:J,3:K+1) - 2 * U(2:J,2:K) + U(2:J,1:K-1) );
   t = t + dt;    count = count + 1;
   fprintf('.')
end
fprintf('\n  completed N = %d steps, average dt = %.4f\n',count,tf/count)

surf(x,y,U),  shading('interp'),  xlabel x,  ylabel y,  zlabel u

