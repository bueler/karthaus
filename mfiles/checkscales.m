% checkscales.m calculates scales d, epsilon, [u], [tau] from
%   measured scales l, [A], [a] and constants rho, g, and n=3,
%   following page 331 of Fowler (1997)
% note: we have already chosen scale [w] := eps [u] for vertical
%   velocity scale, and then notice [w] = [a]

ell = 3000e3;
g = 9.81;
secpera = 31556926.0;

% alternate defns of ice density; first gives "rho g = 0.1 bar m^-1"
%rho = (0.1 * 1e5) / 9.81;
rho = 910;

abracket = 0.1 / secpera;  % 0.1 m a-1

Abracket = 0.2 * (1e5)^(-3) / secpera;  % 0.2 bar^-3 a^-1

d = ( (ell^4 * abracket) / (2 * Abracket * (rho * g)^3) )^(1/8);

epsilon = d / ell;

taubracket = rho * g * d * epsilon;

ubracket = 2 * Abracket * d * taubracket^3;

fprintf('choice of scales gives\n')
fprintf('   d       = %8.3f m\n   epsilon = %8.5f\n',...
        d, epsilon)
fprintf('   [u]     = %8.3f m a-1\n   [tau]   = %8.2f Pa\n',...
        ubracket * secpera, taubracket)

