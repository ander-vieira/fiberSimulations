function [lightPout, electricPout] = dyeInCladRaytracing(dopant, N, diameter, q, lightL, darkL, incidenceAngleDegrees)
%DYEINCLADRAYTRACING Simulate fibers using raytracing as the main tool
%   This function simulates a fiber to obtain a resulting power output by
%   running "photons" through the fiber using raytracing, as opposed to
%   using rate equations.
%   This function simulates two interphase fibers with the dye dopant being
%   in the fiber's CLADDING only (as opposed to dyeRaytracing).
%
%   dopant: the dye dopant's name (see getDyeDopantAttributes)
%   N: the dopant concentration (molecules/m^3)
%   diameter: the fiber's diameter (m)
%   q: the fraction of the core's diameter and the fiber's (Din/Dout)
%   lightL: the fiber length under sunlight (m)
%   darkL: the fiber length NOT under sunlight (e.g. connectors) (m)
%   incidenceAngleDegrees: inclination of the incident sunlight (º)

tic;

Mx = 100;
Mz = 100;
Mn = 200;
M = Mx*Mz*Mn; % Number of photons to simulate

rng('shuffle');

dlambda = 5e-9;
[minLambda, maxLambda] = getLambdaRanges(dopant, dlambda);
ll = minLambda:dlambda:maxLambda;
numll = length(ll);

% Rotation due to sunlight incidence angle
incidenceAngle = deg2rad(incidenceAngleDegrees);
rotationMatrix = [1 0 0 ; 0 cos(incidenceAngle) -sin(incidenceAngle) ; 0 sin(incidenceAngle) cos(incidenceAngle)]';

[tauRad, sigmaabsFun, sigmaemiFun, tauNR] = getDyeDopantAttributes(dopant);

% Get probability distributions for sunlight photons and emitted photons
solarDistribution = solarIrradianceSpline(ll); % W/m^3
solarConstant = sum(solarDistribution); % W/m^3
solarDistribution = solarDistribution/solarConstant;
solarConstant = solarConstant*dlambda; % W/m^2
sigmaabsValues = sigmaabsFun(ll);
sigmaemiValues = sigmaemiFun(ll);
emittedDistribution = sigmaemiValues/sum(sigmaemiValues);
alfaPMMAValues = attenuationPMMA(ll);
nPMMA = refractionIndexPMMA(ll);

quantumYield = tauNR/(tauRad+tauNR);

% Initial values
incomingPower = solarConstant*diameter*lightL*cos(incidenceAngle); % W
Pout = zeros(1, numll); % W
totalPhotons = 0;
finalPhotons = 0;

%%%%%%%%%
% Auxiliary functions

% Using a discrete probability distribution, generate random wavelength
function i = generateDistributedLambda(ll, distribution)
    randomNumber = rand();
    totalSum = 0;
    
    for i = 1:length(ll)
        totalSum = totalSum + distribution(i);
        
        if(randomNumber < totalSum)
            break;
        end
    end
end

function medium = getMedium(position)
    distanceToCenter = position(1)^2+position(2)^2;
    
    if(distanceToCenter < (diameter*q)^2/4)
        % Fiber core
        medium = 0;
    elseif(distanceToCenter < diameter^2/4)
        % Fiber cladding
        medium = 1;
    else
        % Air
        medium = 2;
    end
end

function n = getMediumParams(medium, k)
    if medium == 0
        % Fiber core
        n = nPMMA(k);
    elseif medium == 1
        % Fiber cladding
        n = nPMMA(k)+0.2;
    elseif medium == 2
        % Air
        n = 1;
    end
end

function ds = getAbsorptionDistance(alpha)
    if alpha <= 0
        ds = inf;
    else
        ds = -log(rand)/alpha;
    end
end

function R = fresnelR(N1, N2, cos1, cos2)
    R = (((N1*cos1-N2*cos2)/(N1*cos1+N2*cos2))^2+((N1*cos2-N2*cos1)/(N1*cos2+N2*cos1))^2)/2;
end

% Get the point of intersection with an interphase
function [ds, medium] = getPhaseIntersection(position, direction)
    beta = (direction(1)*position(1, 1)+direction(2)*position(1, 2))/(direction(1)^2+direction(2)^2);
    gamma1 = (position(1, 1)^2+position(1, 2)^2-(diameter*q)^2/4)/(direction(1)^2+direction(2)^2);
    gamma2 = (position(1, 1)^2+position(1, 2)^2-diameter^2/4)/(direction(1)^2+direction(2)^2);
    
    dsList = [sqrt(beta^2-gamma2)-beta -sqrt(beta^2-gamma2)-beta sqrt(beta^2-gamma1)-beta -sqrt(beta^2-gamma1)-beta];
    ds = min(dsList(imag(dsList) == 0 & dsList > 0));
    
    if isempty(ds)
        ds = Inf;
    end
    
    if ds == dsList(1)
        medium = 2;
    elseif ds == dsList(2) || ds == dsList(3)
        medium = 1;
    else
        medium = 0;
    end
end

% Get intersection of ray trajectory with plane given by z = zPos
function ds = getEndIntersection(position, direction, zPos)
    if direction(3) == 0
        % Ray is perpendicular to plane, meaning no intersection
        ds = Inf;
    else
        % Get plane intersection point
        ds = (zPos-position(3))/direction(3);
        
        % Negative intersections don't count
        if ds <= 0
            ds = Inf;
        end
        
        % Check if intersection is inside the fiber diameter
        newPosition = position+direction*ds;
        if newPosition(1)^2+newPosition(2)^2 >= (diameter/2)^2
            ds = Inf;
        end
    end
end

% Run an individual photon through the fiber
% Calls itself recursively for reflected photons and stimulated emission
function runPhoton(position, direction, k, photonPower)
    totalPhotons = totalPhotons + 1;
    
    alphaPMMA = alfaPMMAValues(k);
    alphaDopant = N*sigmaabsValues(k);
    
    while true
        medium = getMedium(position);

        % Get intersection with the fiber's interphases
        [dsPhase, newMedium] = getPhaseIntersection(position, direction);
        
        % Get random distance before ray is absorbed by PMMA
        if medium == 0 || medium == 1
            dsPMMA = getAbsorptionDistance(alphaPMMA);
        else
            dsPMMA = Inf;
        end
        
        % Get random distance before ray is absorbed by dopant
        % Only if ray is currently in cladding
        if medium == 1
            dsDopant = getAbsorptionDistance(alphaDopant);
        else
            dsDopant = Inf;
        end
        
        % Get distance before ray reaches edges of fiber
        dsLeft = getEndIntersection(position, direction, 0);
        dsRight = getEndIntersection(position, direction, lightL+darkL);
        
        % Get distance to the event that happens first
        ds = min([dsPhase dsPMMA dsDopant dsLeft dsRight]);
        
        if ds == Inf
            % Ray escaped the fiber
            break;
        end
        
        position = position+ds*direction;
        
        if ds == dsPMMA
            % Ray was absorbed by PMMA and was lost
            break;
        elseif ds == dsLeft
            % Ray left fiber through left end and was lost
            
            % Reflect ray by a planar mirror on the left end
            direction(3) = -direction(3);
            
            % break;
        elseif ds == dsRight
            % Ray left fiber through right end and was concentrated
            
            % Add power carried by the ray to total output power
            Pout(k) = Pout(k) + photonPower;
            finalPhotons = finalPhotons + 1;
            
            break;
        elseif ds == dsDopant
            % Ray was absorbed by dopant
            
            if rand > quantumYield
                % Ray didn't get reemitted and was lost
                break;
            else
                % Ray was reemitted via spontaneous emission
                
                % Generate random direction
                rand1 = 2*pi*rand();
                rand2 = 2*rand()-1;
                direction = [sqrt(1-rand2^2)*cos(rand1) sqrt(1-rand2^2)*sin(rand1) rand2];
                
                % Generate new wavelength after emission
                oldLambda = ll(k);
                k = generateDistributedLambda(ll, emittedDistribution);
                alphaPMMA = alfaPMMAValues(k);
                alphaDopant = N*sigmaabsValues(k);
                
                % Change photon power (Stokes shift loss)
                photonPower = photonPower*oldLambda/ll(k);
            end
        elseif ds == dsPhase
            % Ray reached interphase

            nIndex = getMediumParams(medium, k);
            nIndexNew = getMediumParams(newMedium, k);
            
            % Get incidence angle
            normalVector = [position(1) position(2) 0];
            normalVector = normalVector/vecnorm(normalVector);

            cosTheta = normalVector*direction';
            sinTheta = sqrt(1-cosTheta^2);
            
            % Apply Snell's law
            newSinTheta = sinTheta*nIndex/nIndexNew;

            reflectedDirection = direction-2*cosTheta*normalVector;

            if newSinTheta >= 1
                % Total internal reflection (TIR)
                direction = reflectedDirection;
            else
                % Refraction or partial reflection
                newCosTheta = sign(cosTheta)*sqrt(1-newSinTheta^2);

                R = fresnelR(nIndex, nIndexNew, cosTheta, newCosTheta);

                if rand < R
                    % Partial reflection
                    direction = reflectedDirection;
                else
                    % Refraction
                    
                    % Calculate direction after refraction vis Snell's law
                    tangentVector = direction-cosTheta*normalVector;
                    newTangentVector = tangentVector*nIndex/nIndexNew;
                    newNormalVector = newCosTheta*normalVector;
                    refractedDirection = newTangentVector+newNormalVector;
                    
                    direction = refractedDirection;
                end
            end

            % Slightly displace the ray off the interphase
            position = position + direction*1e-9;
        end
    end
end

%%%%%%%%%
% Main loop

if nargout == 0
    fprintf('Progress: 000%%');
end

index = 0;
% Get direction of ray (rotated using incidence angle)
direction = [0 -1 0]*rotationMatrix;

for i = 1:Mx
    posX = diameter*((i-1/2)/Mx-1/2);
    
    for j = 1:Mz
        posZ = lightL*(j-1/2)/Mz;

        % Get initial position of ray (directly above fiber, random)
        position = [posX 0 posZ]+[0 (0.6*diameter)/rotationMatrix(2, 2) 0]*rotationMatrix;
        
        for n = 1:Mn
            % Generate random photon from incident sunlight
            k = generateDistributedLambda(ll, solarDistribution);
            photonPower = incomingPower/M;

            runPhoton(position, direction, k, photonPower);

            index = index + 1;

            if nargout == 0 && mod(index*100, M) == 0
                fprintf('\b\b\b\b%03d%%', index*100/M);
            end
        end
    end
end

if nargout == 0
    fprintf('\n');
end

%%%%%%%%%
% Display results

lightPout = sum(Pout);

diodeSurface = pi*diameter^2/4; % m^2
electricPout = solarCellConversion(ll, Pout, diodeSurface);

if nargout == 0
    fprintf('Simulation time: %.1f s\n', toc());
    fprintf('Output power of fiber: %g uW\n', lightPout*1e6);
    fprintf('Output power of solar cell: %g uW\n', electricPout*1e6);
    fprintf('Photons reaching output: %d/%d\n', finalPhotons, totalPhotons);
    
    % Plot power spectrum
    figure(1);
    plot(ll*1e9, Pout*1e-3/dlambda);
    title('Power spectrum at end of fiber (raytracing method)');
    xlabel('\lambda (nm)');
    ylabel('Power spectrum (\muW/nm)');
end

end