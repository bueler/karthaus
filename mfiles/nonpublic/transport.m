function transport(J,tfyears,method,problem)
% TRANSPORT   Artificial problems to test transport schemes.  Equation is
%   H_t = M - (u H)_x
% on  0 <= x <= L,  where L = 100 km, to compute H(t,x).  Uses boundary value
% H(t,0) = 2000 m.  The velocity is fixed and time-independent:  u(x) = 100 (x/L)^2.
% form:
%   transport(J,tfyears,method,problem)
% where:
%   J = (number of grid subintervals)
%   tfyears = (final time, in years)
%   method = (1 for upwind, 2 for Lax-Wendroff)
%   problem = (1 for constant mass balance  M = 0.5 m/a, 
%              2 for mass balance from manufactured solution)
% example to compare methods:
%   >> transport(1000,1000,1,1)
%   >> transport(1000,1000,2,1)
% example to evaluate errors made by Lax-Wendroff:
%   >> transport(100,100,2,2)
%   >> transport(300,100,2,2)
%   >> transport(1000,100,2,2)
%   >> transport(3000,100,2,2)

if (method==1)
  fprintf('  TRANSPORT  using upwind method ')
else
  fprintf('  TRANSPORT  using Lax-Wendroff method ')
end
if (problem==1)
  fprintf('on problem with spatially-constant mass balance\n')
else
  fprintf('on problem with manufactured mass balance\n')
end

secpera = 3.1556926e7;
L = 1.0e5; % m
M0 = 0.5 / secpera;  % constant mass balance for problem 1
Hleft = 2000.0;
tf = tfyears * secpera;

dx = L/J;
x = 0:dx:L; % length J+1

H0 = Hleft - 1.0e-2 * x;

u_fcn = @(x) 100.0 * (x/L).^2 / secpera;
u = u_fcn(x);
if (problem == 2)
  ux_fcn = @(x) (200.0 / L) * (x/L) / secpera;
end
%figure(2), plot(x/1000,u * secpera), xlabel('x  (km)'), ylabel('velocity  (m/a)')
if (method == 2)
  xinter = (dx/2):dx:L-(dx/2);  % staggered grid locations; length J
  uinter = u_fcn(xinter);  % time-independent staggered grid velocities
end

dt = 0.5 * dx / max(u);  % enforce CFL-ish:   max(u) * dt / dx = 1/2
N = ceil(tf / dt);  dt = tf / N;
fprintf('  TRANSPORT  doing %d steps of %.4f years ...\n',N,dt/secpera)

M = M0 * ones(size(x));  % setting mass balance if problem == 1
if (method == 2)
  Minter = M0 * ones(size(xinter)); 
  Mhalf = M0 * ones(size(x));
end

H = H0;
figure(1), plot(x/1000,H), hold on
maxH = max(H);
for n=1:N
  if (problem == 2)
    t = n * dt;
    M = Mmanufactured(x,t,L,Hleft,secpera,u_fcn,ux_fcn);
    if (method == 2)
      Minter = Mmanufactured(xinter,t,L,Hleft,secpera,u_fcn,ux_fcn);
      Mhalf = Mmanufactured(x,t+(dt/2),L,Hleft,secpera,u_fcn,ux_fcn);
    end
  end
  if (method == 1)  % upwind
    H(2:J+1) = H(2:J+1) + dt * (M(2:J+1) - ( u(2:J+1) .* H(2:J+1) - u(1:J) .* H(1:J) ) / dx);
  else % Lax-Wendroff
    Hinter = 0.5 * (H(1:J) + H(2:J+1));
    Hinter = Hinter + (dt/2) * (Minter - ( u(2:J+1) .* H(2:J+1) - u(1:J) .* H(1:J) ) / dx);
    H(J+1) = H(J+1) + dt * (M(J+1) - ( u(J+1) * H(J+1) - u(J) * H(J) ) / dx);
    H(2:J) = H(2:J) + dt * (Mhalf(2:J) - ...
         ( uinter(2:J) .* Hinter(2:J) - uinter(1:J-1) .* Hinter(1:J-1) ) / dx);
  end
  maxH = max(maxH,max(H));
end
plot(x/1000,H,'r')
hold off, grid on, xlabel('x  (km)'), ylabel('thickness  (m)')
fprintf('  TRANSPORT  done.  result:  H(t_f,L) = %.3f m\n',H(end))

if (problem == 2)
  err = H - Hmanufactured(x,t,L,Hleft,secpera);
  fprintf('  TRANSPORT  reports max error = %.6f m\n',max(abs(err)))
end

  function [H, Ht, Hx] = Hmanufactured(x,t,L,Hleft,secpera)
    H0 = Hleft - 1.0e-2 * x;
    H = H0 - 300 * (x/L) * (t / (100*secpera));
    Ht = - 300 * (x/L) / (100*secpera);
    Hx = -1.0e-2 - (300 / L) * (t / (100*secpera));

  function M = Mmanufactured(x,t,L,Hleft,secpera,u_fcn,ux_fcn)
    u = u_fcn(x);
    ux = ux_fcn(x);
    [H, Ht, Hx] = Hmanufactured(x,t,L,Hleft,secpera);
    M = Ht + ux .* H + u .* Hx;

