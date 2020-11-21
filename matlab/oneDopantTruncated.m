function [lightPout, electricPout] = oneDopantTruncated(dopant, N, diameter, lightL, darkL)
%ONEDOPANTITERATIVE Summary of this function goes here
%   This function is cursed. Please don't use it.

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

[tau, sigmaabsFun, sigmaemiFun] = getDyeDopantAttributes(dopant);
sigmaabs = sigmaabsFun(ll);
sigmaemi = sigmaemiFun(ll);
wnsp = sigmaemi / sum(sigmaemi);
truncatedWnsp = wnsp.*(ll>ll');
truncatedWnsp = truncatedWnsp./(sum(truncatedWnsp, 2)+realmin);

alfaPMMA = valuesalfaPMMA(ll);

isol = solarIrradianceSpline(ll);

ncore = refractionIndexPMMA(ll);
beta = (ncore - 1)./(2*ncore);

efficiency = zeros(1, numll);
for k = 1:numll
    alfaCore = alfaPMMA(k) + sigmaabs(k)*N + realmin;
%     efficiency(k) = fiberAbsorptionNoReflections(ncore(k), diameter, sigmaabs(k)*N, alfaCore);
%     efficiency(k) = fiberAbsorptionReflections(ncore(k), diameter, sigmaabs(k)*N, alfaCore);
    efficiency(k) = fiberAbsorptionTwoInterfaces(ncore(k), 1.4, diameter, .98, sigmaabs(k)*N, alfaCore, alfaPMMA(k));
end

% Precalculated constants
concentrationToPower = pi*h*c*diameter^2./(4*ll);
Nsolconst = diameter*isol*dlambda.*efficiency./concentrationToPower;
Nabsconst = sigmaabs./concentrationToPower;
Nestconst = sigmaemi./concentrationToPower;
Pattconst = (alfaPMMA+N*sigmaabs)*dz;
PNconst1 = concentrationToPower.*beta.*truncatedWnsp*dz/tau;
PNconst2 = (sigmaabs+sigmaemi)*dz;

P = zeros(numzz, numll);
Pleft = zeros(numzz, numll);
N2 = zeros(numzz, numll);

error = 1;
while error > 1e-8
    % Boundary condition for P
    P(1, :) = zeros(1, numll);
    Pleft(end, :) = zeros(1, numll);
    
    previousP = P(end, :);
    
    % Update N2
    for j = 1:numzz-1
        for k = 1:numll
            evalP = (P(j, :)+P(j+1, :)+Pleft(j, :)+Pleft(j+1, :))/2;
            wabs = Nabsconst(k)*evalP(k);
            west = sum(Nestconst.*evalP.*[zeros(1, k) ones(1, numll-k)]);

            A = 1/tau+wabs+west;

            if j <= lightj
                b = Nsolconst(k)+N*wabs;
            else
                b = N*wabs;
            end

            N2(j, k) = A\b;
        end
    end
    
    % Update P
    for j = 1:numzz-1
        for k = 1:numll
            evalP = P(j, k);
            
            P(j+1, k) = P(j, k);
            
            P(j+1, k) = P(j+1, k) - Pattconst(k)*evalP;
            
            for m = 1:numll
                evalN2 = N2(j, m);
                
                P(j+1, k) = P(j+1, k) + PNconst1(m, k)*evalN2;
                P(j+1, k) = P(j+1, k) + PNconst2(k)*evalN2*evalP;
            end
        end
    end
    
    % Update Pleft
    for jinv = numzz:-1:2
        for k = 1:numll
            evalP = Pleft(jinv, k);
            
            Pleft(jinv-1, k) = Pleft(jinv, k);
            
            Pleft(jinv-1, k) = Pleft(jinv-1, k) - Pattconst(k)*evalP;
            
            for m=1:numll
                evalN2 = N2(jinv-1, m);
                
                Pleft(jinv-1, k) = Pleft(jinv-1, k) + PNconst1(m, k)*evalN2;
                Pleft(jinv-1, k) = Pleft(jinv-1, k) + PNconst2(k)*evalN2*evalP;
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

if nargout == 0
    fprintf('Simulation time: %.1f s\n', toc());
    fprintf('Output power of fiber: %g uW\n', lightPout*1e6);
    fprintf('Output power of solar cell: %g uW\n', electricPout*1e6);
        
    figure(1);
    plot(ll*1e9, Pout*1e-3/dlambda);
    title('Power spectrum at end of fiber (FDM method)');
    xlabel('\lambda (nm)');
    ylabel('Power spectrum (\muW/nm)');
end

end

