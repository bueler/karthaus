function H = upwind(x,ustag,H0,Hleft,M,tf,firstorder)
% UPWIND  Solves
%    H_t + (uH)_x = M
% by upwinding, both second-order (default) and first order.
% Assumes H(x=0) = Hleft given.  Assumes u(x=0) is positive and
% that u(x(end)) is positive, Dirichlet condition on
% H is only need at left (x=0) end.
% form:
%   FIXME
% where ...
% Example:  TESTUPWIND

spera = 31556926.0;
J = length(x);  dx = x(2)-x(1);
dtCFL = dx / max(abs(ustag));
N = ceil(tf/dtCFL);  dt = tf/N;
fprintf('dt for CFL = %.3f a;  doing N = %d steps to final time tf = %.3f a\n',
        dtCFL/spera,N,tf/spera);

H = H0;
for k = 1:N
  Hnew(2:J) = H(2:J) + dt * (q(   FIXME
  H(1) = Hleft;
end
