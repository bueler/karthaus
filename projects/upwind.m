function H = upwind(x,u,Hinitial,H0,M,tf)
% UPWIND  Solves
%    H_t + (u H)_x = M
% by conservative first order upwinding (i.e. the donor cell method).
% Assumes H(x=0) = H0 given.  Assumes u(x=0) is positive, and that
% u(x(end)) is positive, so that Dirichlet condition on H is only needed
% at left (x=0) end.
% Form:
%   H = upwind(x,u,Hinitial,H0,M,tf)
% where
%   H        = output; ice thickness at tf (s) run; J+1 vector (m)
%   x        = grid; must be equally spaced; J+1 vector (m)
%   u        = velocity is assumed time-independent; J+1 vector (m/s)
%   Hinitial = initial ice thickness; J+1 vector (m)
%   H0       = left hand ice thickness boundary condition, at x(1); scalar (m)
%   M        = surface mass balance; J+1 vector (m/s)
%   tf       = duration of run; scalar (s)
% Note x,u,Hinitial,M,H are all on the "regular" grid even though
% internally the velocity and flux are computed on the "staggered" grid.
% Uses the CFL condition to determine the time step.  As text output,
% one '.' is printed per time step.
% Example:  TESTUPWIND

spera = 31556926.0;
J = length(x)-1;                 % J+1 points  x_1,x_2,...,x_{J+1}
dx = x(2)-x(1);
dtCFL = dx / max(abs(u));
N = ceil(tf/dtCFL);
dt = tf/N;  % divide into integer # of steps
nu = dt / dx;
fprintf('dt for CFL = %.3f a;  doing N = %d steps to final time tf = %.3f a\n',
        dtCFL/spera,N,tf/spera);

ustag = (u(1:J) + u(2:J+1))/2;
qstag = zeros(size(ustag));  % merely pre-allocate
H = Hinitial;
for k = 1:N
  % compute fluxes at cell boundaries
  for j = 1:J
    if ustag(j) >= 0
      qstag(j) = ustag(j) * H(j);
    else
      qstag(j) = ustag(j) * H(j+1);
    end
  end
  % make time-step: update H
  Hnew(1)   = H0;
  Hnew(2:J) = H(2:J) + dt * M(2:J) - nu * (qstag(2:J) - qstag(1:J-1));
  Hnew(J+1) = H(J+1) + dt * M(J+1) - 2 * nu * (u(J+1) * H(J+1) - qstag(J));
  H = Hnew;  % actually update
  fprintf('.')
end
fprintf('\n')
