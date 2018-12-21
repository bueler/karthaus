% SHELFCONV  show convergence study of SSA solver
% calls: TESTSHELF

J = [50 100 200 400 800 1600 3200];

dxkm = 200.0 ./ J;
maxerr = ones(size(dxkm));  % allocate
for j=1:length(J)
  fprintf('J = %4d:\n',J(j))
  [av,maxerr(j)] = testshelf(J(j));
  maxerr(j) = 3.1556926e7 * maxerr(j);
  fprintf('  max error = %.5f m/a\n',maxerr(j))
end

figure(1)  % show convergence plot
loglog(dxkm,maxerr,'o-','markersize',16,'linewidth',2.0)
pf = polyfit(log(dxkm),log(maxerr),1);
hold on, loglog(dxkm,exp(pf(1)*log(dxkm)+pf(2)),'r:','linewidth',2.0), hold off
grid on, xlabel('dx  (km)','fontsize',20)
ylabel('maximum error  (m/a)','fontsize',20)
result = sprintf('rate O(dx^{%.2f})\n',pf(1));
disp(result)
text(2*dxkm(end-1),maxerr(end-1),result,'fontsize',20,'color','r')

%figure(2)  % on super-low resolution grid, show solution
%testshelf(20);

