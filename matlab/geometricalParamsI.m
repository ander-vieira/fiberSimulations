function [beta, Kz] = geometricalParamsI(ncore)
%GEOMETRICALPARAMSI Summary of this function goes here
%   Detailed explanation goes here

m = sqrt(1-1/ncore^2);

bracket1 = @(a) asin(m./sin(a))-(m./sin(a)).*sqrt(1-m^2./sin(a).^2);
integrand1 = @(a) sin(a).*bracket1(a);
integrand2 = @(a) sin(a).*cos(a).*bracket1(a);

beta = ((ncore-1)/ncore + 2/pi*quadgk(integrand1, asin(m), pi/2))/2;
Kz = (2*beta)/(m^2/2+2/pi*quadgk(integrand2, asin(m), pi/2));

end

