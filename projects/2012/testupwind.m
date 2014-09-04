function maxdiff = testupwind(J)
% TESTUPWIND  Tests the UPWIND code by setting up the van der Veen
% exact ice shelf solution and then treating it as time dependent.
% Any changes are a measure of the numerical error in the scheme for
% the mass continuity equation, i.e. in the scheme used in UPWIND.
% This code is a place to start for idea 1 in ideas.txt.
% Examples:
%   >> testupwind           % uses J = 100 subintervals; about 5 sec
%   >> testupwind(50)       % try coarser
%   >> testupwind(500)      % try finer; takes a few minutes
% Calls:  UPWIND, EXACTSHELF

if nargin < 1, J = 100; end

L = 100e3;        % m
M0pera = 0.3;     % m/a
Hg = 500.0;       % m
ugpera = 50.0;    % m/a

x = linspace(0,L,J+1);
spera = 31556926.0;
[uexact,Hinitial] = exactshelf(x,L,M0pera/spera,Hg,ugpera/spera); % exact steady
M = (M0pera/spera) * ones(size(x));

% short run
Hshort  = upwind(x,uexact,Hinitial,Hg,M,  100.0*spera);
% long run
Hlong   = upwind(x,uexact,Hinitial,Hg,M,  300.0*spera);
% longer run
Hlonger = upwind(x,uexact,Hinitial,Hg,M,30000.0*spera);

figure(1), clf
subplot(4,1,[1 2])
plot(x/1000.0,Hinitial,x/1000.0,Hshort,x/1000.0,Hlong,x/1000.0,Hlonger)
legend('Hinitial','H at 100 a','H at 300 a','H at 30000 a')
xlabel('x  (km)'),  ylabel('thickness  (m)')
subplot(4,1,3)
plot(x/1000.0,(Hshort-Hinitial),x/1000.0,(Hlong-Hinitial),...
     x/1000.0,(Hlonger-Hinitial))
legend('at 100 a','at 300 a','at 30000 a')
xlabel('x  (km)'),  ylabel('thickness difference from initial  (m)')
subplot(4,1,4),  plot(x/1000.0,uexact*spera)
xlabel('x  (km)'),  ylabel('velocity  (m/a)')

maxdiff = max([max(abs(Hshort-Hinitial)) max(abs(Hlong-Hinitial)) max(abs(Hlonger-Hinitial))]);
