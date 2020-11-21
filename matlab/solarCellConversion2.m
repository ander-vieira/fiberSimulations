function [electricPout] = solarCellConversion2(ll, lightPin, diodeSurface)

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

responsivity = photodiodeResponsivity(ll);

diodeCurrent = sum(lightPin .* responsivity);

% a and b are defined to simplify the equation to solve
b2 = Isat2 / Isat1;
b = diodeCurrent / Isat1;

% Equation to solve: 1 + b2 + bL = (2a+1)*exp(2a) + b2*(a+1)*exp(a)

a = realmin; % Avoids division by 0
while 1
    expMinusA = exp(-a);
    deltaA = ((1+b2+b)*expMinusA^2 - (2*a+1) - b2*(a+1)*expMinusA)/(2*(1+b2+b)*expMinusA^2+2-b2*a*expMinusA);
%     deltaA = ((1+b2+b)*expMinusA^2-(2*a+1)-b2*(a+1)*expMinusA)/(4*(a+1)+b2*(a+2)*expMinusA);
    a = a + deltaA;
    
    if abs(deltaA/a)<1e-6
        break;
    end
end

Vmax = 2*thermalVoltage*a;
Imax = diodeCurrent - Isat1*(exp(2*a)-1) - Isat2*(exp(a));
electricPout = Vmax*Imax;

end

