% TESTUPWIND  does what you think

spera = 31556926.0;
J = 101;  L = 100e3;
x = linspace(0,L,J);
M0 = 0.3/spera;
Hg = 500;  ug = 50/spera;
[uexact,H0] = exactshelf(x,L,M0,Hg,ug); % exact steady
u = uexact;  M = M0*ones(size(x));
H = upwind(x,u,H0,M,10.0*spera);

plot(x/1000.0,H), xlabel('x  (km)'),  ylabel('thickness  (m)')
