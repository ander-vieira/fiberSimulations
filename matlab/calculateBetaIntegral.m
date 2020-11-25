function beta = calculateBetaIntegral(lambda)
%CALCULATEBETAINTEGRAL Summary of this function goes here
%   Detailed explanation goes here

ncore = refractionIndexPMMA(lambda);
m = sqrt(1-1/ncore^2);

integrand = @(a) (asin(sqrt(1-m^2./sin(a).^2))+m./sin(a).*sqrt(1-m^2./sin(a).^2)).*sin(a);

% integrand = @(a) asin(sqrt(1-m^2./sin(a).^2))+m./sin(a).*sqrt(1-m^2./sin(a).^2);

beta = (1 - 2/pi*quadgk(integrand, asin(m), pi/2))/2;

end

