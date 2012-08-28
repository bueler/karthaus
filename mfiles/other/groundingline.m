FIXME:  needs to get result from buildant.m

% compute grounding line from Antarctic 50km gridded data
%   set, derived from SeaRISE-Antarctica data, then plot 
%   ice profile (bed elevation and surface elevation) along it
% preparation:
%   >> help buildant
%   >> buildant
% uses fields x,y,thk,topg

rhoi = 910.0;
rhow = 1028.0;
grounded = rhoi * thk + rhow * topg;  % if this is positive then
                                      %   rho_i H > - rho_w b
                                      % so ice is grounded

% show "raw" contour:
figure(97)
C = contour(x/1000 , y/1000, grounded, [0, 0]);  % specify zero contour
axis equal, xlabel('x (km)'), ylabel('y (km)')

% now inspect C, which is 2 x N, to get polygonal grounding line
% see "help contourc" for format
% note contour C has several closed loops, and we grab the *biggest* one;
%   the rest are islands and lakes, etc.
jpoly = 1;
Npoly = C(2,1);
j = 1;
while j + C(2,j) < size(C,2)
  j = j + C(2,j) + 1;
  if Npoly < C(2,j)
    jpoly = j;
    Npoly = C(2,j);
  end
end

% we know what part of C is the main contour
xpoly = C(1,jpoly+1:jpoly+Npoly);  % note units of xpoly,ypoly are km
ypoly = C(2,jpoly+1:jpoly+Npoly);

figure(98)
plot(xpoly, ypoly)
axis equal, xlabel('x (km)'), ylabel('y (km)')

% make polygonal grounding line M times as fine before interpolating;
%   xpf = "x polygonal fine", so to speak
M = 10;  
finelength = (length(xpoly)-1) * M + 1;
xpf = zeros(1,finelength);
ypf = zeros(1,finelength);
for j=1:length(xpoly)-1
  xpf((j-1)*M+1:j*M) = linspace(xpoly(j),xpoly(j+1),M);
  ypf((j-1)*M+1:j*M) = linspace(ypoly(j),ypoly(j+1),M);
end
xpf(end) = xpoly(end);
ypf(end) = ypoly(end);
figure(99), plot(xpf, ypf), axis square

% compute arclength sf (= s fine) along grounding line
sf = zeros(1,finelength);
sf(2:finelength) = cumsum(sqrt(diff(xpf).^2 + diff(ypf).^2));
L = max(sf);
disp(['computed length of grounding line is L = ' num2str(L) ' km'])

% interpolate to get profile; see "help interp2"
thkf = interp2(x/1000,y/1000,thk,xpf,ypf,'linear');

% plot profile
figure(99)
r = rhoi/rhow;
plot(sf, (1-r) * thkf, sf, -r * thkf)
xlabel('arclength (km)'), ylabel('profile in (m)')
xtloc = linspace(0,L,6);
set(gca,'XTick',xtloc)
set(gca,'XTickLabel',['0    '; '0.2 L'; '0.4 L'; '0.6 L'; '0.8 L'; 'L    '])

% profile detail along Siple coast and Ross area
figure(100)
sfe = sf(sf>0.8*L);
thkfe = thkf(sf>0.8*L);
plot(sfe, (1-r) * thkfe, sfe, -r * thkfe)
xlabel('arclength (km)'), ylabel('profile in (m)')
set(gca,'XTick',[0.8*L, L])
set(gca,'XTickLabel',['0.8 L'; 'L    '])

% go back to map-plane curve and label locations
figure(98)
hold on
text(xpf(1),ypf(1),'s = 0 = L')
for k=1:4
  i = max(find(sf< (k/5)*L));
  text(xpf(i),ypf(i),['s = ' num2str(k/5) 'L'])
end
hold off

