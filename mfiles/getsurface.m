function h = getsurface(H,b,rho,rhow)
% GETSURFACE  helper function to build surface elevations carefully,
% according to grounded or floating
% example:  see usage inside siageneral.m

f = rho / rhow;                     % fraction of floating ice below surface
h = H + b;                          % only valid where H > 0 and grounded
h(H <= 0) = max(b(H <= 0),0.0);     % if no ice
floating = (H > 0) & (b < - f * H); % points where flotation criterion
h(floating) = (1-f) * H(floating);  %   ... is applied

