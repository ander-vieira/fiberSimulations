function [lightPout, electricPout] = earthDopantFdm(dopant, N, diameter, lightL, darkL)
%EARTHDOPANTFDM Summary of this function goes here
%   Detailed explanation goes here

tic;

if nargin < 5
    darkL = 0;
end

% Constants
c = 3e8; % Speed of light (m/s)
h = 6.63e-34; % Planck constant (J*s)

minlambda = 340e-9;
dlambda = 2e-9;
maxlambda = 740e-9;

dt = 3e-7;
dz = 5e-5; % m

zz = 0:dz:(lightL+darkL);
numzz = length(zz);
lightj = lightL/dz;

ll = minlambda:dlambda:maxlambda;
numll = length(ll);

[tauT, tauD, wTD, wDT, sigmaabsFun, sigmaemiFun] = getEarthDopantAttributes(dopant);
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
NTespconst = dt*wDT/(wDT+wTD)/tauT;
NDespconst = dt*wTD/(wDT+wTD)/tauD;
Nsolconst = sum(isol*dlambda*dt*diameter.*efficiency./concentrationToPower);
Nabsconst = Kz.*sigmaabs*dt./concentrationToPower;
Nestconst = Kz.*sigmaemi*wTD/(wDT+wTD)*dt./concentrationToPower;
Ppropconst = ncore*dz/(c*dt);
Pattconst = Kz.*alfaPMMA*dz;
Pabsconst1 = Kz.*sigmaabs*dz;
Pabsconst2 = Kz.*wDT/(wDT+wTD)*sigmaabs*dz;
Pestconst = Kz.*sigmaemi*wTD/(wDT+wTD)*dz;
Pespconst = concentrationToPower.*beta.*wnsp*dz*wTD/(wDT+wTD)/tauD;

P = zeros(2, numzz, numll);
Pleft = zeros(2, numzz, numll);
NZ = zeros(2, numzz-1);

i = 0;
imin = 100;
error = 0;
while i < imin || error > 1e-7
    P(1, :, :) = P(2, :, :);
    Pleft(1, :, :) = Pleft(2, :, :);
    NZ(1, :) = NZ(2, :);
    
    % Boundary condition for P
    P(2, 1, :) = zeros(1, 1, numll);
    Pleft(2, end, :) = zeros(1, 1, numll);
    
    % Update NZ (NT and ND together)
    for j = 1:numzz-1
        evalNZ = NZ(1, j);
        
        NZ(2, j) = NZ(1, j);
        
        NZ(2, j) = NZ(2, j) - evalNZ*NTespconst;
        NZ(2, j) = NZ(2, j) - evalNZ*NDespconst;
        
        if j <= lightj
            NZ(2, j) = NZ(2, j) + Nsolconst;
        end
        
        for k = 1:numll
            evalP = (P(1, j, k)+P(1, j+1, k)+Pleft(1, j, k)+Pleft(1, j+1, k))/2;
            
            NZ(2, j) = NZ(2, j) + Nabsconst(k)*evalP*(N-evalNZ);
            NZ(2, j) = NZ(2, j) - Nestconst(k)*evalP*evalNZ;
        end
    end
    
    % Update P
    for j = 1:numzz-1
        evalNZ = NZ(2, j);
        
        for k = 1:numll
            evalP = P(2, j, k);
            evaldP = P(2, j, k)-P(1, j, k);
            
            P(2, j+1, k) = P(2, j, k);
            
            P(2, j+1, k) = P(2, j+1, k) - Ppropconst(k)*evaldP;
            P(2, j+1, k) = P(2, j+1, k) - Pattconst(k)*evalP;
            
            P(2, j+1, k) = P(2, j+1, k) + Pespconst(k)*evalNZ;
            P(2, j+1, k) = P(2, j+1, k) + Pestconst(k)*evalP*evalNZ;
            P(2, j+1, k) = P(2, j+1, k) - Pabsconst1(k)*evalP*N;
            P(2, j+1, k) = P(2, j+1, k) + Pabsconst2(k)*evalP*evalNZ;
        end
    end
    
    % Update Pleft
    for jinv = numzz:-1:2
        evalNZ = NZ(2, jinv-1);
        
        for k = 1:numll
            evalP = Pleft(2, jinv, k);
            evaldP = Pleft(2, jinv, k)-Pleft(1, jinv, k);
            
            Pleft(2, jinv-1, k) = Pleft(2, jinv, k);
            
            Pleft(2, jinv-1, k) = Pleft(2, jinv-1, k) - Ppropconst(k)*evaldP;
            Pleft(2, jinv-1, k) = Pleft(2, jinv-1, k) - Pattconst(k)*evalP;
            
            Pleft(2, jinv-1, k) = Pleft(2, jinv-1, k) + Pespconst(k)*evalNZ;
            Pleft(2, jinv-1, k) = Pleft(2, jinv-1, k) + Pestconst(k)*evalP*evalNZ;
            Pleft(2, jinv-1, k) = Pleft(2, jinv-1, k) - Pabsconst1(k)*evalP*(N-evalNZ);
        end
    end
    
    i = i + 1;
    % Calculate the error (for the loop condition)
    error = max(abs((P(2, end, :)-P(1, end, :))./(P(2, end, :)+realmin)));
end

elapsedTime = i*dt;

Pout = squeeze(P(2, end, :))';
lightPout = sum(Pout);

diodeSurface = pi*diameter^2/4; % m^2
electricPout = solarCellConversion(ll, Pout, diodeSurface);

estimatedError = max(NZ(2, :))*diameter*max(sigmaabs);

if nargout == 0
    fprintf('Simulation time: %.1f s\n', toc());
    fprintf('Elapsed time: %.1f us\n', elapsedTime*1e6);
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

