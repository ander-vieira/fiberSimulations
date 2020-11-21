function absorption = fiberAbsorptionNoReflections(ncore, diameter, alfaDopant, alfaCore)
%FIBERABSORPTIONNOREFLECTIONS Absorbed fraction model without reflections
%   Calculates the fraction of incident sunlight absorbed by a fiber
%   without taking any reflections into account.
%   Exact solution via integrating over u

% Distance traveled inside the fiber
distanceInFiber = @(u) diameter*sqrt(1 - u.^2/ncore^2);

absExp = @(u) exp(-alfaCore*distanceInFiber(u));

% Add the fraction that is absorbed by the dopant itself
integrand = @(u) (1 - absExp(u))*alfaDopant/alfaCore;

% Integrate numerically
absorption = quadgk(integrand, 0, 1);

end