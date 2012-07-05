% SHELFCONV  show convergence study of SSA solver using testshelf.m

J = [50 100 200 400 800 1600 3200];
dxkm = 200.0 ./ J;
fprintf('(showing one dot per outer iteration)\n')
for j=1:length(J)
  fprintf('J = %4d,  dx = %5.3f km:\n',J(j),dxkm(j))
  [av,maxerr(j)] = testshelf(J(j));
  maxerr(j) = 3.1556926e7 * maxerr(j);
  fprintf('  max error = %.5f m/a\n',maxerr(j))
end

figure(1)  % show convergence plot
loglog(dxkm,maxerr,'o-','markersize',16)
pf = polyfit(log(dxkm),log(maxerr),1);
hold on, loglog(dxkm,exp(pf(1)*log(dxkm)+pf(2)),'r:'), hold off
grid on, xlabel('dx  (km)','fontsize',16)
ylabel('maximum error  (m a^{-1})','fontsize',16)
result = sprintf('convergence rate O(dx^{%.5f})\n',pf(1));
disp(result),   text(2*dxkm(end-1),maxerr(end-1),result,'fontsize',18)

figure(2)  % on super-low resolution grid, show solution
testshelf(20);  fprintf('\n')

