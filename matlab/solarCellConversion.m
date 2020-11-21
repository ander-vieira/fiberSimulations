function [electricPout] = solarCellConversion(ll, lightPin, diodeSurface)

if nargin == 0
    dlambda = 2e-9;
    minlambda = 300e-9;
    maxlambda = 1100e-9;
    ll = minlambda:dlambda:maxlambda;
    
    diameter = 1e-3;
    diodeSurface = pi*diameter^2/4;
    
    isol = solarIrradianceSpline(ll);
    lightPin = diodeSurface*isol*dlambda;
end

kB = 1.3806e-23; % Boltzmann constant (J/K)
electronCharge = 1.6e-19; % C

Isat1 = 1e-8*diodeSurface;
Isat2 = .5e-3*diodeSurface;
Tamb = 300; % Ambient temperature (K)
thermalVoltage = kB*Tamb/electronCharge;
fillFactor = 0.8;

responsivity = photodiodeResponsivity(ll);

diodeCurrent = sum(lightPin .* responsivity);
openCircuitVoltage = 2*thermalVoltage*log((sqrt((2*Isat1+Isat2)^2+4*Isat1*diodeCurrent)-Isat2)/(2*Isat1));

electricPout = diodeCurrent*openCircuitVoltage*fillFactor;

end

