function [lightPout, electricPout] = combinedDopantIterative(dyeDopant, dyeN, earthDopant, earthN, solarType, diameter, lightL, darkL)
%COMBINEDDOPANTITERATIVE Simulate a POF with dye and rare earth dopants
%   Allows multiple dye dopants AND multiple rare earth dopants
%   Iterative method: gets the stationary state of the fiber
%   but still allows simulating bidirectional propagation
%   
%   Parameters:
%   dyeDopant: array with the names of each dye dopant
%   dyeN: array with the concentration of each dye dopant (m^-3)
%   Must be the same length as dyeDopant
%   earthDopant: array with the names of each rare earth dopant
%   earthN: array with the concentration of each rare earth dopant (m^-3)
%   Must be the same length as earthDopant
%   diameter: the fiber's diameter (m)
%   lightL: illuminated length of the fiber (m)
%   darkL: non-illuminated length at the end of the fiber (m)
%   
%   Example calls:
%   combinedDopantIterative(["C1" "C6"], [.7 1]*5.9226e22, [""], [0], 1e-3, .06, .03);
%   combinedDopantIterative([""], [0], ["AC46"], [3e23], 1e-3, .1);
%   combinedDopantIterative([""], [0], ["AC46" "AC46"], [2e23 1e23], 1e-3, .1);
%   combinedDopantIterative(["LumogenRed"], [1.77e22], ["AC46"], [2.96e22], 1e-3, 0.045, .03);

tic;

if nargin < 8
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

numDyeDopants = length(dyeDopant);
dyeTauRad = zeros(numDyeDopants, 1);
dyeTauNR = zeros(numDyeDopants, 1);
dyeSigmaabs = zeros(numDyeDopants, numll);
dyeSigmaemi = zeros(numDyeDopants, numll);
dyeWnsp = zeros(numDyeDopants, numll);

for m = 1:numDyeDopants
    [dyeTauRad(m), sigmaabsFun, sigmaemiFun, dyeTauNR(m)] = getDyeDopantAttributes(dyeDopant(m));
    
    dyeSigmaabs(m, :) = sigmaabsFun(ll);
    dyeSigmaemi(m, :) = sigmaemiFun(ll);
    dyeWnsp(m, :) = dyeSigmaemi(m, :)/sum(dyeSigmaemi(m, :));
end

numEarthDopants = length(earthDopant);
earthTauT = zeros(numEarthDopants, 1);
earthTauD = zeros(numEarthDopants, 1);
earthwTD = zeros(numEarthDopants, 1);
earthwDT = zeros(numEarthDopants, 1);
earthSigmaabs = zeros(numEarthDopants, numll);
earthSigmaemi = zeros(numEarthDopants, numll);
earthWnsp = zeros(numEarthDopants, numll);

for m = 1:numEarthDopants
    [earthTauT(m), earthTauD(m), earthwTD(m), earthwDT(m), sigmaabsFun, sigmaemiFun] = getEarthDopantAttributes(earthDopant);
    
    earthSigmaabs(m, :) = sigmaabsFun(ll);
    earthSigmaemi(m, :) = sigmaemiFun(ll);
    earthWnsp(m, :) = earthSigmaemi(m, :)/sum(earthSigmaemi(m, :));
end

alfaPMMA = attenuationPMMA(ll);

isol = solarIrradianceSpline(ll, solarType);

ncore = refractionIndexPMMA(ll);

beta = zeros(1, numll);
Kz = zeros(1, numll);
for k = 1:numll
%     [beta(k), Kz(k)] = geometricalParamsB(ncore(k));
    [beta(k), Kz(k)] = geometricalParamsI(ncore(k));
end

dyeEfficiency = zeros(numDyeDopants, numll);
for m = 1:numDyeDopants
    for k = 1:numll
        alfaCore = alfaPMMA(k) + sum(dyeSigmaabs(:, k).*dyeN')+ sum(earthSigmaabs(:, k).*earthN') + realmin;
        dyeEfficiency(m, k) = fiberAbsorptionTwoInterfaces(ncore(k), 1.4, diameter, .98, dyeSigmaabs(m, k)*dyeN(m), alfaCore, alfaPMMA(k));
    end
end

earthEfficiency = zeros(numEarthDopants, numll);
for m = 1:numEarthDopants
    for k = 1:numll
        alfaCore = alfaPMMA(k) + sum(dyeSigmaabs(:, k).*dyeN')+ sum(earthSigmaabs(:, k).*earthN') + realmin;
        earthEfficiency(m, k) = fiberAbsorptionTwoInterfaces(ncore(k), 1.4, diameter, .98, earthSigmaabs(m, k)*earthN(m), alfaCore, alfaPMMA(k));
    end
end

% Precalculated constants
concentrationToPower = pi*h*c*diameter^2./(4*ll);
dyeNsolconst = diameter*sum(isol*dlambda.*dyeEfficiency./concentrationToPower, 2);
earthNsolconst = diameter*sum(isol*dlambda.*earthEfficiency./concentrationToPower, 2);
dyeNabsconst = Kz.*dyeSigmaabs./concentrationToPower;
dyeNestconst = Kz.*dyeSigmaemi./concentrationToPower;
earthNabsconst = Kz.*earthSigmaabs./concentrationToPower;
earthNestconst = Kz.*earthSigmaemi./concentrationToPower;
Pattconst = Kz.*(alfaPMMA+dyeN*dyeSigmaabs+earthN*earthSigmaabs)*dz;
dyePNconst1 = concentrationToPower.*beta.*dyeWnsp*dz./dyeTauRad;
dyePNconst2 = Kz.*(dyeSigmaabs+dyeSigmaemi)*dz;
earthPNconst1 = concentrationToPower.*beta.*earthWnsp*dz./earthTauD;
earthPNconst2 = Kz.*earthSigmaemi*dz;
earthPNconst3 = Kz.*earthSigmaabs*dz;

P = zeros(numzz, numll);
Pleft = zeros(numzz, numll);
N2 = zeros(numzz, numDyeDopants);
NT = zeros(numzz, numEarthDopants);
ND = zeros(numzz, numEarthDopants);

error = 1;
while error > 1e-8
    % Boundary condition for P
    P(1, :) = zeros(1, numll);
    Pleft(end, :) = zeros(1, numll);
    
    % Keep this for calculating the error
    previousP = P(end, :);
    
    for j = 1:numzz-1
        evalP = (P(j, :)+P(j+1, :)+Pleft(j, :)+Pleft(j+1, :))/2;
        
        % Update N2 (dye dopants)
        for m = 1:numDyeDopants
            wabs = sum(dyeNabsconst(m, :).*evalP);
            west = sum(dyeNestconst(m, :).*evalP);
            
            A = 1/dyeTauRad(m)+1/dyeTauNR(m)+wabs+west;
            
            if j <= lightj
                b = dyeNsolconst(m)+dyeN(m)*wabs;
            else
                b = dyeN(m)*wabs;
            end
            
            N2(j, m) = A\b;
        end
        
        % Update NT and ND (rare earth dopants)
        for m = 1:numEarthDopants
            wabs = sum(earthNabsconst(m, :).*evalP);
            west = sum(earthNestconst(m, :).*evalP);

            A = [1/earthTauT(m)+earthwTD(m)+wabs -earthwDT(m) ; -earthwTD(m) 1/earthTauD(m)+earthwDT(m)+west];

            if j <= lightj
                b = [earthNsolconst(m)+earthN(m)*wabs ; 0];
            else
                b = [earthN(m)*wabs ; 0];
            end

            c = A\b;

            NT(j, m) = c(1);
            ND(j, m) = c(2);
        end
    end
    
    % Update P (use FDM over the positive z dimension)
    for j = 1:numzz-1
        for k = 1:numll
            evalP = P(j, k);
            
            P(j+1, k) = P(j, k);
            
            % Attenuation
            P(j+1, k) = P(j+1, k) - Pattconst(k)*evalP;
            
            % Effects of each dye dopant
            for m = 1:numDyeDopants
                evalN2 = N2(j, m);
        
                P(j+1, k) = P(j+1, k) + dyePNconst1(m, k)*evalN2;
                P(j+1, k) = P(j+1, k) + dyePNconst2(m, k)*evalN2*evalP;
            end
            
            % Effects of each earth dopant
            for m = 1:numEarthDopants
                evalNT = NT(j, m);
                evalND = ND(j, m);
            
                P(j+1, k) = P(j+1, k) + earthPNconst1(m, k)*evalND;
                P(j+1, k) = P(j+1, k) + earthPNconst2(m, k)*evalND*evalP;
                P(j+1, k) = P(j+1, k) + earthPNconst3(m, k)*evalNT*evalP;
            end
        end
    end
    
    % Update Pleft (use FDM over the negative z dimension)
    for jinv = numzz:-1:2
        for k = 1:numll
            evalP = Pleft(jinv, k);
            
            Pleft(jinv-1, k) = Pleft(jinv, k);
            
            % Attenuation
            Pleft(jinv-1, k) = Pleft(jinv-1, k) - Pattconst(k)*evalP;
            
            % Effects of each dye dopant
            for m = 1:numDyeDopants
                evalN2 = N2(jinv-1, m);
        
                Pleft(jinv-1, k) = Pleft(jinv-1, k) + dyePNconst1(m, k)*evalN2;
                Pleft(jinv-1, k) = Pleft(jinv-1, k) + dyePNconst2(m, k)*evalN2*evalP;
            end
            
            % Effects of each earth dopant
            for m = 1:numEarthDopants
                evalNT = NT(jinv-1, m);
                evalND = ND(jinv-1, m);
            
                Pleft(jinv-1, k) = Pleft(jinv-1, k) + earthPNconst1(m, k)*evalND;
                Pleft(jinv-1, k) = Pleft(jinv-1, k) + earthPNconst2(m, k)*evalND*evalP;
                Pleft(jinv-1, k) = Pleft(jinv-1, k) + earthPNconst3(m, k)*evalNT*evalP;
            end
        end
    end
    
    % Calculate the error (for the loop condition)
    error = max(abs((P(end, :)-previousP)./(P(end, :)+realmin)));
end

Pout = P(end, :);
lightPout = sum(Pout);

diodeSurface = pi*diameter^2/4; % m^2
electricPout = solarCellConversion2(ll, Pout, diodeSurface);

estimatedError = (sum(max(N2, [], 1)'.*max(dyeSigmaabs,[], 2))+sum(max(NT+ND, [], 1)'.*max(earthSigmaabs, [], 2)))*diameter;

if nargout == 0
    fprintf('Simulation time: %.1f s\n', toc());
    fprintf('Output power of fiber: %g uW\n', lightPout*1e6);
    fprintf('Output power of solar cell: %g uW\n', electricPout*1e6);
    fprintf('Estimated error of approximation: %g\n', estimatedError);
        
    figure(1);
    plot(ll*1e9, Pout*1e-3/dlambda);
    title('Power spectrum at end of fiber (iterative method)');
    xlabel('\lambda (nm)');
    ylabel('Power spectrum (\muW/nm)');
end

end

