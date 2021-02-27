function [electricPout] = solarCellConversion(ll, lightPin, diodeSurface)
%SOLARCELLCONVERSION Obtain electrical power from incoming optical power
%    Uses the entire spectrum of light to apply responsivity
%    Then uses a two-diode model to get V and I, but uses FF to get P

if nargin == 0
    % Use direct sunlight and 1mm circular cell to calculate
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

Isat1 = 1e-8*diodeSurface; % A
Isat2 = .5e-3*diodeSurface; % A
Tamb = 300; % Ambient temperature (K)
thermalVoltage = kB*Tamb/electronCharge; % V
fillFactor = 0.8;

solarCellR = solarCellResponsivity(ll);

diodeCurrent = sum(lightPin .* solarCellR); % I
openCircuitVoltage = 2*thermalVoltage*log((sqrt((2*Isat1+Isat2)^2+4*Isat1*diodeCurrent)-Isat2)/(2*Isat1)); % V

electricPout = diodeCurrent*openCircuitVoltage*fillFactor; % W

end

