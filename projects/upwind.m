function H = upwind(x,u,H0,M,tf,firstorder)
% UPWIND  Solves
%    H_t + (uH)_x = M
% by upwinding, both second-order (default) and first order.
% form:
%   FIXME
% where ...
% Example:  TESTUPWIND

spera = 31556926.0;
J = length(x);  dx = x(2)-x(1);
dtCFL = dx / max(abs(u));
N = ceil(tf/dtCFL);  dt = tf/N;
fprintf('dt for CFL = %.3f a;  doing N = %d steps to final time tf = %.3f a\n',
        dtCFL/spera,N,tf/spera);

H = H0;
for k = 1:N
  % FIXME
end  
