function [lightPout, electricPout, ll, Pout] = dyeFdm(dopant, N, diameter, q, lightL, darkL)
%DYEFDM Summary of this function goes here
%   Detailed explanation goes here

tic;

if nargin < 5
    darkL = 0;
end

% Constants
c = 3e8; % Speed of light (m/s)
h = 6.63e-34; % Planck constant (J*s)

dlambda = 2e-9;
[minLambda, maxLambda] = getLambdaRanges(dopant, dlambda);

dt = 3e-11;
dz = 5e-5; % m

zz = 0:dz:(lightL+darkL);
numzz = length(zz);
lightj = lightL/dz;

ll = minLambda:dlambda:maxLambda;
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
    efficiency(k) = fiberAbsorptionTwoInterfaces(ncore(k), 1.4, diameter, q, sigmaabs(k)*N, alfaCore, alfaPMMA(k));
end

% Precalculated constants
concentrationToPower = pi*h*c*diameter^2./(4*ll);
Nespconst = dt/tauRad+dt/tauNR;
Nsolconst = sum(isol*dlambda*dt*diameter.*efficiency./concentrationToPower);
Nabsconst = Kz.*sigmaabs*dt./concentrationToPower;
Nestconst = Kz.*sigmaemi*dt./concentrationToPower;
Ppropconst = ncore*dz/(c*dt);
Pattconst = Kz.*alfaPMMA*dz;
Pabsconst = Kz.*sigmaabs*dz;
Pestconst = Kz.*sigmaemi*dz;
Pespconst = concentrationToPower.*beta.*wnsp*dz/tauRad;

P = zeros(2, numzz, numll);
Pleft = zeros(2, numzz, numll);
N2 = zeros(2, numzz-1);

% Samples along the time axis for plot over time
sampleDT = 0.3e-9;
numSamples = 10;
curSamples = 1;
Psamples = zeros(numSamples, 1);

i = 0;
imin = 100;
error = 0;
while i < imin || error > 1e-8
    P(1, :, :) = P(2, :, :);
    Pleft(1, :, :) = Pleft(2, :, :);
    N2(1, :) = N2(2, :);
    
    % Boundary condition for P
    P(2, 1, :) = zeros(1, 1, numll);
    Pleft(2, end, :) = zeros(1, 1, numll);
    
    % Update N2
    for j = 1:numzz-1
        evalN2 = N2(1, j);
        
        N2(2, j) = N2(1, j);
        
        N2(2, j) = N2(2, j) - evalN2*Nespconst;
        
        if j <= lightj
            N2(2, j) = N2(2, j) + Nsolconst;
        end
        
        for k = 1:numll
            evalP = (P(1, j, k)+P(1, j+1, k)+Pleft(1, j, k)+Pleft(1, j+1, k))/2;
            
            N2(2, j) = N2(2, j) + Nabsconst(k)*evalP*(N-evalN2);
            N2(2, j) = N2(2, j) - Nestconst(k)*evalP*evalN2;
        end
    end
    
    % Update P
    for j = 1:numzz-1
        evalN2 = N2(2, j);
        
        for k = 1:numll
            evalP = P(2, j, k);
            evaldP = P(2, j, k)-P(1, j, k);
            
            P(2, j+1, k) = P(2, j, k);
            
            P(2, j+1, k) = P(2, j+1, k) - Ppropconst(k)*evaldP;
            P(2, j+1, k) = P(2, j+1, k) - Pattconst(k)*evalP;
            
            P(2, j+1, k) = P(2, j+1, k) + Pespconst(k)*evalN2;
            P(2, j+1, k) = P(2, j+1, k) + Pestconst(k)*evalP*evalN2;
            P(2, j+1, k) = P(2, j+1, k) - Pabsconst(k)*evalP*(N-evalN2);
        end
    end
    
    % Update Pleft
    for jinv = numzz:-1:2
        evalN2 = N2(2, jinv-1);
        
        for k = 1:numll
            evalP = Pleft(2, jinv, k);
            evaldP = Pleft(2, jinv, k)-Pleft(1, jinv, k);
            
            Pleft(2, jinv-1, k) = Pleft(2, jinv, k);
            
            Pleft(2, jinv-1, k) = Pleft(2, jinv-1, k) - Ppropconst(k)*evaldP;
            Pleft(2, jinv-1, k) = Pleft(2, jinv-1, k) - Pattconst(k)*evalP;
            
            Pleft(2, jinv-1, k) = Pleft(2, jinv-1, k) + Pespconst(k)*evalN2;
            Pleft(2, jinv-1, k) = Pleft(2, jinv-1, k) + Pestconst(k)*evalP*evalN2;
            Pleft(2, jinv-1, k) = Pleft(2, jinv-1, k) - Pabsconst(k)*evalP*(N-evalN2);
        end
    end
    
    % Take the time axis samples (dynamic array allocation)
    if(mod(i, floor(sampleDT/dt)) == 0)
        curSamples = curSamples + 1;
        
        if(curSamples > numSamples)
            % Expand array to double the current size
            Psamples = [Psamples ; zeros(numSamples, 1)];
            numSamples = numSamples*2;
        end
        
        % Sample the total output power at end of fiber
        Psamples(curSamples) = sum(P(2, end, :));
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

estimatedError = max(N2(2, :))*diameter*max(sigmaabs);

if nargout == 0
    fprintf('Simulation time: %.1f s\n', toc());
    fprintf('Elapsed time: %.1f ns\n', elapsedTime*1e9);
    fprintf('Output power of fiber: %g uW\n', lightPout*1e6);
    fprintf('Output power of solar cell: %g uW\n', electricPout*1e6);
    fprintf('Estimated error of approximation: %g\n', estimatedError);
        
    % Plot power spectrum
    figure(1);
    plot(ll*1e9, Pout*1e-3/dlambda);
    title('Power spectrum at end of fiber (FDM method)');
    xlabel('\lambda (nm)');
    ylabel('Power spectrum (\muW/nm)');
    
    % Plot over time axis
    figure(2);
    sampleTT = sampleDT*(0:(curSamples-1));
    plot(sampleTT*1e9, Psamples(1:curSamples)*1e6);
end

end

