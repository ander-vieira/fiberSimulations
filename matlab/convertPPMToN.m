function Nout = convertPPMToN(Nin, Mw, rho, reverse)
%CONVERTPPMTON Convert ppm to m^-3 concentrations and viceversa
%   Nin: concentration in ppm/m^-3
%   Mw: molecular weight of dopant (g/mol)
%   rho: density of medium (PMMA: 1.18e6) (g/m^3)
%   reverse: if true, convert m^-3 to ppm instead
%   Nout: concentration in m^-3/ppm

N_A = 6.022e23; % Avogadro's number

if nargin < 3
    rho = 1.18e6; % Default value: density of PMMA
end

if nargin >= 4 && reverse
    Nout = (Nin*Mw/rho/N_A)*1e6; % m^-3 to ppm
else
    Nout = (Nin/1e6)*N_A*rho/Mw; % ppm to m^-3
end

end

