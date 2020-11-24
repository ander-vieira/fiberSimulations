function [electricPout] = solarCellConversion2(ll, lightPin, diodeSurface)
%SOLARCELLCONVERSION Obtain electrical power from incoming optical power
%    Uses the entire spectrum of light to apply responsivity
%    Then uses a two-diode model to get V and I
%    And uses numerical minmax algorithm (Newton-Raphson) to get P

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

responsivity = photodiodeResponsivity(ll);

% Get Iph from the incoming sunlight
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

% a = V/V_T/2
Vmax = 2*thermalVoltage*a; % V
Imax = diodeCurrent - Isat1*(exp(2*a)-1) - Isat2*(exp(a)); % A
electricPout = Vmax*Imax; % W
% Coupled resistance for maximum power: Vmax/Imax (ohm)

end

