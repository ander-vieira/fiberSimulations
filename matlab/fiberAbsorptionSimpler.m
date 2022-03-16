function absorption = fiberAbsorptionSimpler(ncore, nclad, diameter, q, alfaDopant, alfaCore, alfaClad)
%FIBERABSORPTIONREFLECTIONS Absorbed fraction model without reflections
%   Calculates the fraction of incident sunlight absorbed by a fiber
%   with a one interface model for reflections (without cladding)
%   Exact solution via integrating over u

uCut = min(1, q*nclad);

% Distance traveled inside the fiber
dcore = @(u) diameter*sqrt(q^2 - u.^2/ncore^2);
dclad = @(u) diameter*(sqrt(1 - u.^2/nclad^2) - sqrt(q^2 - u.^2/nclad^2));

cosAir = @(u) sqrt(1 - u.^2);
cosClad1 = @(u) sqrt(1 - u.^2/nclad^2);
cosClad2 = @(u) sqrt(1 - u.^2/(q*nclad)^2);
cosCore = @(u) sqrt(1 - u.^2/(q*ncore)^2);

% Fresnel reflection coefficients in each interface
R_F1 = @(u) (((cosAir(u) - nclad*cosClad1(u))./(cosAir(u) + nclad*cosClad1(u))).^2 + ((nclad*cosAir(u) - cosClad1(u))./(nclad*cosAir(u) + cosClad1(u))).^2)/2;
R_F2 = @(u) (((nclad*cosClad2(u) - ncore*cosCore(u))./(nclad*cosClad2(u) + ncore*cosCore(u))).^2 + ((ncore*cosClad2(u) - nclad*cosCore(u))./(ncore*cosClad2(u) + nclad*cosCore(u))).^2)/2;
T_F1 = @(u) 1-R_F1(u);
T_F2 = @(u) 1-R_F2(u);

expCore = @(u) exp(-alfaCore*dcore(u));
expClad = @(u) exp(-alfaClad*dclad(u));

D1 = @(u) (1-R_F2(u).*expCore(u)).*(1-R_F1(u).*R_F2(u).*expClad(u).^2)-R_F1(u).*T_F2(u).^2.*expClad(u).^2.*expCore(u);
D2 = @(u) T_F1(u).*T_F2(u).*expClad(u).*(1-expCore(u))*alfaDopant/alfaCore;

% Result of the infinite reflections is obtained analytically (see notes)
% Add the fraction that is absorbed by the dopant itself
integrand = @(u) D2(u)./D1(u);

% Integrate numerically
absorption = quadgk(integrand, 0, uCut);

end