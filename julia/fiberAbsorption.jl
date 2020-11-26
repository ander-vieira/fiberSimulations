module FiberAbsorption

using QuadGK;

noReflections = function(ncore, diameter, alfaDopant, alfaCore)

distanceInFiber(u) = diameter.*sqrt(1 .- u.^2 ./ncore^2);

absExp(u) = exp(-alfaCore*distanceInFiber(u));

integrand(u) = (1 - absExp(u))*alfaDopant/alfaCore;

absorption = quadgk(integrand, 0, 1);

return absorption[1];

end

end
