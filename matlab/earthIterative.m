function [lightPout, electricPout] = earthIterative(dopant, N, diameter, q, lightL, darkL)
%EARTHITERATIVE Simulate fibers through the stationary state equations
%   This function simulates a fiber's behavior under sunlight by solving
%   the rate equations in their stationary state version (i.e. by nulling
%   all derivatives with respect to time).
%   It also uses an iterative approach to run power in both directions
%   through the fiber.
%
%   dopant: the rare earth dopant's name (see getDyeDopantAttributes)
%   N: the dopant concentration (molecules/m^3)
%   diameter: the fiber's diameter (m)
%   q: the fraction of the core's diameter and the fiber's (Din/Dout)
%   lightL: the fiber length under sunlight (m)
%   darkL: the fiber length NOT under sunlight (e.g. connectors) (m)

tic;

if nargin < 5
    darkL = 0;
end

% Constants
c = 3e8; % Speed of light (m/s)
h = 6.63e-34; % Planck constant (J*s)

dlambda = 2e-9;
[minLambda, maxLambda] = getLambdaRanges(dopant, dlambda);

dz = 5e-5; % m

zz = 0:dz:(lightL+darkL);
numzz = length(zz);
lightj = lightL/dz;

ll = minLambda:dlambda:maxLambda;
numll = length(ll);

[tauT, tauD, wTD, wDT, sigmaabsFun, sigmaemiFun] = getEarthDopantAttributes(dopant);
sigmaabs = sigmaabsFun(ll);
sigmaemi = sigmaemiFun(ll);
wnsp = sigmaemi / sum(sigmaemi);

alfaPMMA = attenuationPMMA(ll);

isol = solarIrradianceSpline(ll);

nPMMA = refractionIndexPMMA(ll);

beta = zeros(1, numll);
Kz = zeros(1, numll);
for k = 1:numll
%     [beta(k), Kz(k)] = geometricalParamsB(ncore(k));
    [beta(k), Kz(k)] = geometricalParamsI(nPMMA(k));
end

efficiency = zeros(1, numll);
for k = 1:numll
    alfaCore = alfaPMMA(k) + sigmaabs(k)*N + realmin;
%     efficiency(k) = fiberAbsorptionNoReflections(nPMMA(k), diameter, sigmaabs(k)*N, alfaCore);
%     efficiency(k) = fiberAbsorptionReflections(nPMMA(k), diameter, sigmaabs(k)*N, alfaCore);
    efficiency(k) = fiberAbsorptionTwoInterfaces(nPMMA(k), 1.4, diameter, q, sigmaabs(k)*N, alfaCore, alfaPMMA(k));
end

% Precalculated constants
concentrationToPower = pi*h*c*diameter^2./(4*ll);
Nsolconst = diameter*sum(isol*dlambda.*efficiency./concentrationToPower);
Nabsconst = Kz.*sigmaabs./concentrationToPower;
Nestconst = Kz.*sigmaemi./concentrationToPower;
Pattconst = Kz.*(alfaPMMA+N*sigmaabs)*dz;
PNconst1 = concentrationToPower.*beta.*wnsp*dz/tauD;
PNconst2 = Kz.*(sigmaabs+sigmaemi)*dz;
PNconst3 = Kz.*sigmaabs*dz;

P = zeros(numzz, numll);
Pleft = zeros(numzz, numll);
NT = zeros(numzz, 1);
ND = zeros(numzz, 1);

error = 1;
while error > 1e-9
    % Boundary condition for P
    P(1, :) = zeros(1, numll);
    Pleft(end, :) = zeros(1, numll);
    
    previousP = P(end, :);
    
    % Update NT and ND
    for j = 1:numzz-1
        evalP = (P(j, :)+P(j+1, :)+Pleft(j, :)+Pleft(j+1, :))/2;
        wabs = sum(Nabsconst.*evalP);
        west = sum(Nestconst.*evalP);
        
        A = [1/tauT+wTD+wabs -wDT+wabs ; -wTD 1/tauD+wDT+west];
        
        if j <= lightj
            b = [Nsolconst+N*wabs ; 0];
        else
            b = [N*wabs ; 0];
        end
        
        c = A\b;
        
        NT(j) = c(1);
        ND(j) = c(2);
    end
    
    % Update P
    for j = 1:numzz-1
        evalNT = NT(j);
        evalND = ND(j);
        
        for k = 1:numll
            evalP = P(j, k);
            
            P(j+1, k) = P(j, k);
            
            P(j+1, k) = P(j+1, k) - Pattconst(k)*evalP;
            
            P(j+1, k) = P(j+1, k) + PNconst1(k)*evalND;
            P(j+1, k) = P(j+1, k) + PNconst2(k)*evalND*evalP;
            P(j+1, k) = P(j+1, k) + PNconst3(k)*evalNT*evalP;
        end
    end
    
    % Update Pleft
    for jinv = numzz:-1:2
        evalNT = NT(jinv-1);
        evalND = ND(jinv-1);
        
        for k = 1:numll
            evalP = Pleft(jinv, k);
            
            Pleft(jinv-1, k) = Pleft(jinv, k);
            
            Pleft(jinv-1, k) = Pleft(jinv-1, k) - Pattconst(k)*evalP;
            
            Pleft(jinv-1, k) = Pleft(jinv-1, k) + PNconst1(k)*evalND;
            Pleft(jinv-1, k) = Pleft(jinv-1, k) + PNconst2(k)*evalND*evalP;
            Pleft(jinv-1, k) = Pleft(jinv-1, k) + PNconst3(k)*evalNT*evalP;
        end
    end
    
    % Calculate the error (for the loop condition)
    error = max(abs((P(end, :)-previousP)./(P(end, :)+realmin)));
end

Pout = P(end, :);
lightPout = sum(Pout);

diodeSurface = pi*diameter^2/4; % m^2
electricPout = solarCellConversion(ll, Pout, diodeSurface);

estimatedError = max(NT+ND)*diameter*max(sigmaabs);

if nargout == 0
    fprintf('Simulation time: %.1f s\n', toc());
    fprintf('Output power of fiber: %g uW\n', lightPout*1e6);
    fprintf('Output power of solar cell: %g uW\n', electricPout*1e6);
    fprintf('Estimated error of approximation: %g\n', estimatedError);
        
    % Plot power spectrum
    figure(1);
    plot(ll*1e9, Pout*1e-3/dlambda);
    title('Power spectrum at end of fiber (iterative method)');
    xlabel('\lambda (nm)');
    ylabel('Power spectrum (\muW/nm)');
end

end

