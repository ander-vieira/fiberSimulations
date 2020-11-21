function [lightPout, electricPout] = multiDopantStationary(dopant, N, diameter, lightL, darkL)
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

numDopants = length(dopant);
tau = zeros(numDopants, 1);
sigmaabs = zeros(numDopants, M);
sigmaemi = zeros(numDopants, M);
wnsp = zeros(numDopants, M);

for m = 1:numDopants
    [tau(m), sigmaabsFun, sigmaemiFun] = getDyeDopantAttributes(dopant(m));
    
    sigmaabs(m, :) = sigmaabsFun(ll);
    sigmaemi(m, :) = sigmaemiFun(ll);
    wnsp(m, :) = sigmaemi(m, :)/sum(sigmaemi(m, :));
end

alfaPMMA = valuesalfaPMMA(ll);

isol = solarIrradianceSpline(ll);

ncore = refractionIndexPMMA(ll);
beta = (ncore - 1)./(2*ncore);

efficiency = zeros(numDopants, M);
for m = 1:numDopants
    for k = 1:M
        alfaCore = alfaPMMA(k) + sum(sigmaabs(:, k).*N') + realmin;
    %     efficiency(m, k) = fiberAbsorptionNoReflections(ncore(k), diameter, sigmaabs(m, k)*N, alfaCore);
    %     efficiency(m, k) = fiberAbsorptionReflections(ncore(k), diameter, sigmaabs(m, k)*N, alfaCore);
        efficiency(m, k) = fiberAbsorptionTwoInterfaces(ncore(k), 1.4, diameter, .98, sigmaabs(m, k)*N(m), alfaCore, alfaPMMA(k));
    end
end

equivalentIrr = sum(isol.*dlambda.*efficiency.*ll, 2);

C = sum(equivalentIrr*diameter*beta.*wnsp./ll, 1);

G = zeros(M);
for i = 1:M
    for j = 1:M
        if(i == j)
            G(i, i) = -sum(sigmaabs(:, i).*N'.*(1-beta(i).*wnsp(:, i)))-alfaPMMA(i)-realmin;
        else
            G(i, j) = sum(sigmaabs(:, j).*N'*beta(i).*wnsp(:, i))*ll(j)/ll(i);
        end
    end
end

Pfinal = expm(G*darkL)*(expm(G*lightL)-eye(M))*(G\C');

N2final = 4*equivalentIrr.*tau/(pi*h*c*diameter)+ 4*N'.*tau.*sum(Pfinal'.*sigmaabs.*ll, 2)/(pi*h*c*diameter^2);
estimatedError = sum(N2final.*max(sigmaabs, [], 2))*diameter;

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

