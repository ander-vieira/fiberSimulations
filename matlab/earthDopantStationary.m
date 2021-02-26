function [lightPout, electricPout] = earthDopantStationary(dopant, N, diameter, lightL, darkL)
%ONEDOPANTSIM_STATIONARY Summary of this function goes here
%   Detailed explanation goes here

tic;

if nargin < 5
    darkL = 0;
end

% Constants
c = 3e8; % Speed of light (m/s)
h = 6.63e-34; % Planck constant (J*s)

minlambda = 240e-9;
dlambda = 2e-9;
maxlambda = 740e-9;

ll = minlambda:dlambda:maxlambda;
M = length(ll);

[tauT, tauD, wTD, wDT, sigmaabsFun, sigmaemiFun] = getEarthDopantAttributes(dopant);
sigmaabs = sigmaabsFun(ll);
sigmaemi = sigmaemiFun(ll);
wnsp = sigmaemi / sum(sigmaemi);

alfaPMMA = valuesalfaPMMA(ll);

isol = solarIrradianceSpline(ll);

ncore = refractionIndexPMMA(ll);

beta = zeros(1, numll);
Kz = zeros(1, numll);
for k = 1:numll
%     [beta(k), Kz(k)] = geometricalParamsB(ncore(k));
    [beta(k), Kz(k)] = geometricalParamsI(ncore(k));
end

qD = wTD/(wTD+wDT+1/tauD);
eqTau  = 1/((1-qD)/tauT+qD/tauD);
R = eqTau*qD/tauD;

efficiency = zeros(1, M);
for k = 1:M
    alfaCore = alfaPMMA(k) + sigmaabs(k)*N + realmin;
%     efficiency(i) = fiberAbsorptionNoReflections(ncore(i), diameter, sigmaabs(i)*N, sigmaabs(i)*N+alfaPMMA(k)+realmin);
%     efficiency(i) = fiberAbsorptionReflections(ncore(i), diameter, sigmaabs(i)*N, alfaPMMA(i));
    efficiency(k) = fiberAbsorptionTwoInterfaces(ncore(k), 1.4, diameter, .98, sigmaabs(k)*N, alfaCore, alfaPMMA(k));
end

equivalentIrr = sum(isol.*dlambda.*efficiency.*ll);

C = equivalentIrr*diameter*beta.*wnsp./ll*R;

G = zeros(M);
for i = 1:M
    for j = 1:M
        if(i == j)
            G(i, i) = -Kz(i)*sigmaabs(i).*N.*(1-beta(i).*wnsp(i)*R)-Kz(i).*alfaPMMA(i)-realmin;
        else
            G(i, j) = Kz(j)*sigmaabs(j)*N*beta(i)*wnsp(i)*R*ll(j)/ll(i);
        end
    end
end

Pfinal = expm(G*darkL)*(expm(G*lightL)-eye(M))*(G\C');

N2final = 4*equivalentIrr*eqTau/(pi*h*c*diameter)+ 4*N*eqTau*sum(Pfinal'.*sigmaabs.*ll)/(pi*h*c*diameter^2);
estimatedError = N2final.*max(sigmaabs)*diameter;

lightPout = sum(Pfinal);
diodeSurface = pi*diameter^2/4; % m^2
electricPout = solarCellConversion(ll, Pfinal', diodeSurface);

if nargout == 0
    fprintf('Simulation time: %.1f s\n', toc());
    fprintf('Output power of fiber: %g uW\n', lightPout*1e6);
    fprintf('Output power of solar cell: %g uW\n', electricPout*1e6);
    fprintf('Estimated error of approximation: %g\n', estimatedError);
    
    figure(1);
    plot(ll*1e9, Pfinal*1e-3/dlambda);
    title('Power spectrum at end of fiber (stationary method)');
    xlabel('\lambda (nm)');
    ylabel('Power spectrum (\muW/nm)');
end

end

