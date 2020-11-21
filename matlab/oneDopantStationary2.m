function [lightPout, electricPout] = oneDopantStationary(dopant, N, diameter, lightL, darkL)
%ONEDOPANTSIM_STATIONARY Summary of this function goes here
%   I have no idea how this works anymore
%   But it's not completely exact when simulating both directions at once
%   May be due to precision since the matrices involved are absurd

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

C = equivalentIrr*diameter*beta.*wnsp./ll;

C = [C -C];

G = zeros(2*M);
for i = 1:2*M
    for j = 1:2*M
        overlapI = mod(i-1, M)+1;
        overlapJ = mod(j-1, M)+1;
        
        if(i == j)
            G(i, i) = -sigmaabs(overlapI).*N.*(1-beta(overlapI).*wnsp(overlapI))-alfaPMMA(overlapI)-realmin;
        else
            G(i, j) = sigmaabs(overlapJ)*N*beta(overlapI)*wnsp(overlapI)*ll(overlapJ)/ll(overlapI);
        end
        
        if i > M
            G(i, j) = -G(i, j);
        end
    end
end

% G = zeros(M);
% for i = 1:M
%     for j = 1:M
%         overlapI = mod(i-1, M)+1;
%         overlapJ = mod(j-1, M)+1;
%         
%         if(i == j)
%             G(i, i) = -sigmaabs(overlapI).*N.*(1-beta(overlapI).*wnsp(overlapI))-alfaPMMA(overlapI)-realmin;
%         else
%             G(i, j) = sigmaabs(overlapJ)*N*beta(overlapI)*wnsp(overlapI)*ll(overlapJ)/ll(overlapI);
%         end
%         
%         if i > M
%             G(i, j) = -G(i, j);
%         end
%     end
% end

[U, D] = eig(G);
v = (U\C')';

d = diag(D)';

T = U;
T(M+1:end, :) = T(M+1:end, :).*exp(d*lightL);

s = sum(U.*v./d, 2);

q = (T\s)';

% Pfinal = sum(U.*(exp(d*lightL).*q-v./d), 2);

Pfinal = zeros(2*M, 1);
for k = 1:2*M
    if d(k) < 0
        Pfinal = Pfinal + U(:, k)*(exp(d(k)*lightL)*q(k)-v(k)/d(k))*exp(d(k)*darkL);
    else
        Pfinal = Pfinal + U(:, k)*(exp(d(k)*lightL)*q(k)-v(k)/d(k))*exp(-d(k)*darkL);
    end
end

Pfinal = Pfinal(1:M);

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

