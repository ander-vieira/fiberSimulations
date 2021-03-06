function [lightPout, electricPout] = multiDyeIterative(dopant, N, diameter, q, lightL, darkL)
%MULTIDYEITERATIVE Summary of this function goes here
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

dz = 5e-5; % m

zz = 0:dz:(lightL+darkL);
numzz = length(zz);
lightj = lightL/dz;

ll = minLambda:dlambda:maxLambda;
numll = length(ll);

numDopants = length(dopant);
tauRad = zeros(numDopants, 1);
tauNR = zeros(numDopants, 1);
sigmaabs = zeros(numDopants, numll);
sigmaemi = zeros(numDopants, numll);
wnsp = zeros(numDopants, numll);

for m = 1:numDopants
    [tauRad(m), sigmaabsFun, sigmaemiFun, tauNR(m)] = getDyeDopantAttributes(dopant(m));
    
    sigmaabs(m, :) = sigmaabsFun(ll);
    sigmaemi(m, :) = sigmaemiFun(ll);
    wnsp(m, :) = sigmaemi(m, :)/sum(sigmaemi(m, :));
end

alfaPMMA = attenuationPMMA(ll);

isol = solarIrradianceSpline(ll);

nPMMA = refractionIndexPMMA(ll);

beta = zeros(1, numll);
Kz = zeros(1, numll);
for k = 1:numll
%     [beta(k), Kz(k)] = geometricalParamsB(ncore(k));
    [beta(k), Kz(k)] = geometricalParamsI(nPMMA(k));
end

efficiency = zeros(numDopants, numll);
for m = 1:numDopants
    for k = 1:numll
        alfaCore = alfaPMMA(k) + sum(sigmaabs(:, k).*N') + realmin;
    %     efficiency(m, k) = fiberAbsorptionNoReflections(ncore(k), diameter, sigmaabs(m, k)*N, alfaCore);
    %     efficiency(m, k) = fiberAbsorptionReflections(ncore(k), diameter, sigmaabs(m, k)*N, alfaCore);
        efficiency(m, k) = fiberAbsorptionTwoInterfaces(nPMMA(k), 1.4, diameter, q, sigmaabs(m, k)*N(m), alfaCore, alfaPMMA(k));
    end
end

% Precalculated constants
concentrationToPower = pi*h*c*diameter^2./(4*ll);
Nsolconst = diameter*sum(isol*dlambda.*efficiency./concentrationToPower, 2);
Nabsconst = Kz.*sigmaabs./concentrationToPower;
Nestconst = Kz.*sigmaemi./concentrationToPower;
Pattconst = Kz.*(alfaPMMA+N*sigmaabs)*dz;
PNconst1 = concentrationToPower.*beta.*wnsp*dz./tauRad;
PNconst2 = Kz.*(sigmaabs+sigmaemi)*dz;

P = zeros(numzz, numll);
Pleft = zeros(numzz, numll);
N2 = zeros(numzz, numDopants);

error = 1;
while error > 1e-8
    % Boundary condition for P
    P(1, :) = zeros(1, numll);
%     P(1, :) = Pleft(1, :); % Mirror
    Pleft(end, :) = zeros(1, numll);
    
    previousP = P(end, :);
    
    % Update N2
    for j = 1:numzz-1
        evalP = (P(j, :)+P(j+1, :)+Pleft(j, :)+Pleft(j+1, :))/2;

        for m = 1:numDopants
            wabs = sum(Nabsconst(m, :).*evalP);
            west = sum(Nestconst(m, :).*evalP);
            
            A = 1/tauRad(m)+1/tauNR(m)+wabs+west;
            
            if j <= lightj
                b = Nsolconst(m)+N(m)*wabs;
            else
                b = N(m)*wabs;
            end
            
            N2(j, m) = A\b;
        end
    end
    
    % Update P
    for j = 1:numzz-1
        for k = 1:numll
            evalP = P(j, k);
            
            P(j+1, k) = P(j, k);
            
            P(j+1, k) = P(j+1, k) - Pattconst(k)*evalP;
            
            for m = 1:numDopants
                evalN2 = N2(j, m);
        
                P(j+1, k) = P(j+1, k) + PNconst1(m, k)*evalN2;
                P(j+1, k) = P(j+1, k) + PNconst2(m, k)*evalN2*evalP;
            end
        end
    end
    
    % Update Pleft
    for jinv = numzz:-1:2
        for k = 1:numll
            evalP = Pleft(jinv, k);
            
            Pleft(jinv-1, k) = Pleft(jinv, k);
            
            Pleft(jinv-1, k) = Pleft(jinv-1, k) - Pattconst(k)*evalP;
            
            
            for m = 1:numDopants
                evalN2 = N2(jinv-1, m);
        
                Pleft(jinv-1, k) = Pleft(jinv-1, k) + PNconst1(m, k)*evalN2;
                Pleft(jinv-1, k) = Pleft(jinv-1, k) + PNconst2(m, k)*evalN2*evalP;
            end
        end
    end
    
    % Calculate the error (for the loop condition)
    error = max(abs((P(end, :)-previousP)./(P(end, :)+realmin)));
end

Pout = P(end, :);
lightPout = sum(Pout);

diodeSurface = pi*diameter^2/4; % m^2
electricPout = solarCellConversion(ll, Pout, diodeSurface);

estimatedError = sum(max(N2, [], 1)'.*max(sigmaabs,[], 2))*diameter;

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

