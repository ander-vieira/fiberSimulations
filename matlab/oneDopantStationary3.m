function [lightPout, electricPout] = oneDopantStationary(dopant, N, diameter, lightL, darkL)
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

[tau, sigmaabsFun, sigmaemiFun] = getDyeDopantAttributes(dopant);
ncore = refractionIndexPMMA(ll);
sigmaabs = sigmaabsFun(ll);
sigmaemi = sigmaemiFun(ll);
wnsp = sigmaemi / sum(sigmaemi);
isol = solarIrradianceSpline(ll);
alfaPMMA = valuesalfaPMMA(ll);

beta = (ncore - 1)./(2*ncore);

efficiency = zeros(1, M);
for k = 1:M
    alfaCore = alfaPMMA(k) + sigmaabs(k)*N + realmin;
%     efficiency(k) = fiberAbsorptionNoReflections(ncore(k), diameter, sigmaabs(k)*N, alfaCore);
%     efficiency(k) = fiberAbsorptionReflections(ncore(k), diameter, sigmaabs(k)*N, alfaCore);
    efficiency(k) = fiberAbsorptionTwoInterfaces(ncore(k), 1.4, diameter, .98, sigmaabs(k)*N, alfaCore, alfaPMMA(k));
end

equivalentIrr = sum(isol.*dlambda.*efficiency.*ll);

C = [zeros(1, M) -2*equivalentIrr*diameter*beta.*wnsp./ll];

G = zeros(2*M);
for i = 1:M
    for j = 1:M
        if(i == j)
            G(i+M, i) = sigmaabs(i).*N.*(1-2*beta(i).*wnsp(i))+alfaPMMA(i)+realmin;
            G(i, i+M) = sigmaabs(i).*N+alfaPMMA(i)+realmin;
        else
            G(i+M, j) = -2*sigmaabs(j)*N*beta(i)*wnsp(i)*ll(j)/ll(i)-realmin;
        end
    end
end

Pinf = -G\C';

Qinf = Pinf(1:M);
Rinf = Pinf(M+1:2*M);
vec = [zeros(M, 1) ; Rinf-Qinf];
mat = [[eye(M) zeros(M)]*(eye(2*M)-expm(G*lightL)) ; eye(M) -eye(M)];

k = mat\vec;

Pdoble = expm(G*darkL)*(Pinf + expm(G*lightL)*k);

Pfinal = (Pdoble(1:M)-Pdoble(M+(1:M)))/2;

N2final = 4*equivalentIrr*tau/(pi*h*c*diameter)+ 4*N*tau*sum(Pfinal'.*sigmaabs.*ll)/(pi*h*c*diameter^2);
estimatedError = N2final*diameter*max(sigmaabs);

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

