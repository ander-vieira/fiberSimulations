function absorption = fiberAbsorptionReflections(ncore, diameter, alfaDopant, alfaCore)
%FIBERABSORPTIONREFLECTIONS Absorbed fraction model without reflections
%   Calculates the fraction of incident sunlight absorbed by a fiber
%   with a one interface model for reflections (without cladding)
%   Exact solution via integrating over u

% Distance traveled inside the fiber
distanceInFiber = @(u) diameter*sqrt(1 - u.^2/ncore^2);

absExp = @(u) exp(-alfaCore*distanceInFiber(u));

cosAir = @(u) sqrt(1-u.^2);
cosCore = @(u) sqrt(1-u.^2/ncore^2);

% Fresnel reflection coefficient in the air-fiber interface
R_F = @(u) (((cosAir(u) - ncore*cosCore(u))./(cosAir(u) + ncore*cosCore(u))).^2+((cosCore(u)-ncore*cosAir(u))./(cosCore(u)+ncore*cosAir(u))).^2)/2;

% Result of the infinite reflections is obtained analytically (see notes)
% Add the fraction that is absorbed by the dopant itself
integrand = @(u) (1 - R_F(u)).*(1 - absExp(u))./(1 - R_F(u) .* absExp(u))*alfaDopant/alfaCore;

% Integrate numerically
absorption = quadgk(integrand, 0, 1);

end