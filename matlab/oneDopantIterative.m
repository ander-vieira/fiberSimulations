function [lightPout, electricPout] = oneDopantIterative(dopant, N, diameter, lightL, darkL)
%ONEDOPANTITERATIVE Summary of this function goes here
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

dz = 5e-5; % m

zz = 0:dz:(lightL+darkL);
numzz = length(zz);
lightj = lightL/dz;

ll = minlambda:dlambda:maxlambda;
numll = length(ll);

[tauRad, sigmaabsFun, sigmaemiFun, tauNR] = getDyeDopantAttributes(dopant);
sigmaabs = sigmaabsFun(ll);
sigmaemi = sigmaemiFun(ll);
wnsp = sigmaemi / sum(sigmaemi);

alfaPMMA = attenuationPMMA(ll);

isol = solarIrradianceSpline(ll);

ncore = refractionIndexPMMA(ll);

beta = zeros(1, numll);
Kz = zeros(1, numll);
for k = 1:numll
%     [beta(k), Kz(k)] = geometricalParamsB(ncore(k));
    [beta(k), Kz(k)] = geometricalParamsI(ncore(k));
end

efficiency = zeros(1, numll);
for k = 1:numll
    alfaCore = alfaPMMA(k) + sigmaabs(k)*N + realmin;
%     efficiency(k) = fiberAbsorptionNoReflections(ncore(k), diameter, sigmaabs(k)*N, alfaCore);
%     efficiency(k) = fiberAbsorptionReflections(ncore(k), diameter, sigmaabs(k)*N, alfaCore);
    efficiency(k) = fiberAbsorptionTwoInterfaces(ncore(k), 1.4, diameter, .98, sigmaabs(k)*N, alfaCore, alfaPMMA(k));
end

% Precalculated constants
concentrationToPower = pi*h*c*diameter^2./(4*ll);
Nsolconst = diameter*sum(isol*dlambda.*efficiency./concentrationToPower);
Nabsconst = Kz.*sigmaabs./concentrationToPower;
Nestconst = Kz.*sigmaemi./concentrationToPower;
Pattconst = Kz.*(alfaPMMA+N*sigmaabs)*dz;
PNconst1 = concentrationToPower.*beta.*wnsp*dz/tauRad;
PNconst2 = Kz.*(sigmaabs+sigmaemi)*dz;

P = zeros(numzz, numll);
Pleft = zeros(numzz, numll);
N2 = zeros(numzz, 1);

error = 1;
while error > 1e-8
    % Boundary condition for P
    P(1, :) = zeros(1, numll);
    Pleft(end, :) = zeros(1, numll);
    
    previousP = P(end, :);
    
    % Update N2
    for j = 1:numzz-1
        evalP = (P(j, :)+P(j+1, :)+Pleft(j, :)+Pleft(j+1, :))/2;
        wabs = sum(Nabsconst.*evalP);
        west = sum(Nestconst.*evalP);
        
        A = 1/tauRad+1/tauNR+wabs+west;
        
        if j <= lightj
            b = Nsolconst+N*wabs;
        else
            b = N*wabs;
        end
        
        N2(j) = A\b;
    end
    
    % Update P
    for j = 1:numzz-1
        evalN2 = N2(j);
        
        for k = 1:numll
            evalP = P(j, k);
            
            P(j+1, k) = P(j, k);
            
            P(j+1, k) = P(j+1, k) - Pattconst(k)*evalP;
            P(j+1, k) = P(j+1, k) + PNconst1(k)*evalN2;
            P(j+1, k) = P(j+1, k) + PNconst2(k)*evalN2*evalP;
        end
    end
    
    % Update Pleft
    for jinv = numzz:-1:2
        evalN2 = N2(jinv-1);
        
        for k = 1:numll
            evalP = Pleft(jinv, k);
            
            Pleft(jinv-1, k) = Pleft(jinv, k);
            
            Pleft(jinv-1, k) = Pleft(jinv-1, k) - Pattconst(k)*evalP;
            Pleft(jinv-1, k) = Pleft(jinv-1, k) + PNconst1(k)*evalN2;
            Pleft(jinv-1, k) = Pleft(jinv-1, k) + PNconst2(k)*evalN2*evalP;
        end
    end
    
    % Calculate the error (for the loop condition)
    error = max(abs((P(end, :)-previousP)./(P(end, :)+realmin)));
end

Pout = P(end, :);
lightPout = sum(Pout);

diodeSurface = pi*diameter^2/4; % m^2
electricPout = solarCellConversion(ll, Pout, diodeSurface);

estimatedError = max(N2)*diameter*max(sigmaabs);

if nargout == 0
    fprintf('Simulation time: %.1f s\n', toc());
    fprintf('Output power of fiber: %g uW\n', lightPout*1e6);
    fprintf('Output power of solar cell: %g uW\n', electricPout*1e6);
    fprintf('Estimated error of approximation: %g\n', estimatedError);
        
    figure(1);
    plot(ll*1e9, Pout*1e-3/dlambda);
    title('Power spectrum at end of fiber (FDM method)');
    xlabel('\lambda (nm)');
    ylabel('Power spectrum (\muW/nm)');
end

end

