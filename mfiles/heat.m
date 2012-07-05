function U = heat(D,J,K,dt,N)
% explicit method for heat equation
%   u_t = D (u_xx + u_yy)
% on -1 < x < 1, -1 < y < 1, for 0 < t < N * dt
% fixed time steps (compare heatadapt.m)
% uses specific gaussian initial condition
% form
%   U = heat(D,J,K,dt,N)
% where
%   D   = diffusivity coeff
%   J,K = number of points in x,y directions
%   dt  = fixed time step
%   N   = number of time steps; final time is tf = N * dt
% try: heat(1.0,30,30,0.001,0);    % initial condition
%      heat(1.0,30,30,0.001,20);   % final time 0.02
%      heat(1.0,30,30,0.004,5);    % final time 0.02; huh?

% spatial grid and initial condition:
dx = 2 / J;    dy = 2 / K;
[x,y] = meshgrid(-1:dx:1, -1:dy:1); % (J+1) x (K+1) grid in x,y plane
U = exp(-30*(x.*x + y.*y));
fprintf('  doing N = %d steps of dt = %.5f for 0.0 < t < %.3f\n',N,dt,N*dt)

% explicit time steps
nu_x = dt * D / (dx*dx);    nu_y = dt * D / (dy*dy);
for n=1:N
   U(2:J,2:K) = U(2:J,2:K) + ...
       nu_x * ( U(3:J+1,2:K) - 2 * U(2:J,2:K) + U(1:J-1,2:K) ) + ...
       nu_y * ( U(2:J,3:K+1) - 2 * U(2:J,2:K) + U(2:J,1:K-1) );
end

surf(x,y,U),  shading('interp'),  xlabel x,  ylabel y,  zlabel u

