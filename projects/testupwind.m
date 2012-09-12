% TESTUPWIND  does what you think

spera = 31556926.0;
J = 101;  L = 100e3;
x = linspace(0,L,J);
M0 = 0.3/spera;
Hg = 500;  ug = 50/spera;
[uexact,H0] = exactshelf(x,L,M0,Hg,ug); % exact steady
u = uexact;  M = M0*ones(size(x));
xstag = (x(1:J-1) + x(2:J)) / 2.0;
ustag = (u(1:J-1) + u(2:J)) / 2.0;
H = upwind(x,ustag,H0,Hg,M,100.0*spera);

subplot(211)
plot(x/1000.0,H),  xlabel('x  (km)'),  ylabel('thickness  (m)')
subplot(212)
plot(xstag/1000.0,ustag*spera),  xlabel('x  (km)'),  ylabel('velocity  (m/a)')
