% using PETSc 3.5.2

% takes about 15 seconds total

% for M in 10 100 1000 10000 100000 1000000; do ./ssaflowline -snes_rtol 1e-12 -ksp_rtol 1e-12 -ssa_epsilon 1e-10 -da_grid_x $M |'grep' 'relative maximum'; done
%                                2.9556e-02 relative maximum
%                                1.0402e-03 relative maximum
%                                1.1305e-05 relative maximum
%                                1.1299e-07 relative maximum
%                                1.1304e-09 relative maximum
%                                4.8689e-11 relative maximum

C = [
10       2.9556e-02;
100      1.0402e-03;
1000     1.1305e-05;
10000    1.1299e-07;
100000   1.1304e-09;
1000000  4.8689e-11];

figure(1)
loglog(C(:,1),C(:,2),'o','markersize',15)
p = polyfit(log(C(2:5,1)),log(C(2:5,2)),1);
fprintf('relative maximum numerical errors go as O(dx^%.3f)\n',-p(1))
hold on
loglog(C(:,1),exp(p(1)*log(C(:,1))+p(2)),'r--')
hold off
legend('error data','fit')
xlabel M, ylabel('relative maximum numerical error'), grid on
print -dpdf numerr.pdf

