function [Pout] = oneDopantRaytracingClad(dopant, N, diameter, q, lightL, darkL, incidenceAngleDegrees)
%ONEDOPANTRAYTRACING2 Simulate fibers using raytracing as the main tool
%   This function simulates a fiber to obtain a resulting power output by
%   running "photons" through the fiber using raytracing, as opposed to
%   using rate equations
%   This function simulates fibers with both core and cladding, with the
%   dopant being only in the core

tic;

M = 1000000; % Number of photons to simulate

rng('shuffle');

% Constants
c = 3e8; % m/s
h = 6.63e-34; % J*s

% Simulation parameters
dz = 1e-4;
zz = 0:dz:(lightL+darkL);
numzz = length(zz);

dlambda = 5e-9;
[minLambda, maxLambda] = getLambdaRanges(dopant, dlambda);
ll = minLambda:dlambda:maxLambda;
numll = length(ll);

da = 1e-5;

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

quantumYield = tauNR/(tauRad+tauNR);
% Conversion constant from power to N2
conversionN2 = 4*ll/(pi*h*c*diameter^2*q^2*dz)/(1/tauRad+1/tauNR); % m^-3/W

% Initial values
N2 = zeros(1, numzz);
incomingPower = solarConstant*diameter*lightL*cos(incidenceAngle); % W
minimumPower = incomingPower/M*0.01;
Pout = zeros(1, numll); % W
totalPhotons = 0;
finalPhotons = 0;
runawayPhotons = 0;
stimulatedPhotons = 0;

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
        n = refractionIndexPMMA(ll(k));
    elseif medium == 1
        % Fiber cladding
        n = 1.4;
    elseif medium == 2
        % Air
        n = 1;
    end
end

% Get the point of intersection with an interphase
function [ds, medium] = getRefractionPoint(position, direction)
    beta = (direction(1)*position(1, 1)+direction(2)*position(1, 2))/(direction(1)^2+direction(2)^2);
    gamma1 = (position(1, 1)^2+position(1, 2)^2-(diameter*q)^2/4)/(direction(1)^2+direction(2)^2);
    gamma2 = (position(1, 1)^2+position(1, 2)^2-diameter^2/4)/(direction(1)^2+direction(2)^2);
    
    dsList = [sqrt(beta^2-gamma2)-beta sqrt(beta^2-gamma1)-beta -sqrt(beta^2-gamma2)-beta -sqrt(beta^2-gamma1)-beta];
    ds = min(dsList(imag(dsList) == 0 & dsList > 0));
    
    if ds == dsList(1)
        medium = 2;
    elseif ds == dsList(2) || ds == dsList(3)
        medium = 1;
    else
        medium = 0;
    end
end

% Run an individual photon through the fiber
% Calls itself recursively for reflected photons and stimulated emission
function runPhoton(position, direction, k, photonPower)
    totalPhotons = totalPhotons + 1;
    
    alfaPMMA = alfaPMMAValues(k);
    sigmaabs = sigmaabsValues(k);
    sigmaemi = sigmaemiValues(k);
    
    position = [zeros(1, 3) ; position];
    
    loopOn = true;
    while loopOn
        position(1, :) = position(2, :);
        
        % Check previous iteration medium
        prevMedium = getMedium(position(1, :));
        prevNIndex = getMediumParams(prevMedium, k);
        
        % Distance interval travelled 
        ds = prevNIndex*da;
        
        % Calculate new position
        position(2, :) = position(1, :) + direction*ds;
        
        % Check current iteration medium
        medium = getMedium(position(2, :));
        
        if prevMedium ~= medium
            % Refraction in an interphase
            
            boundaryDistance = 0.01*da;
            
            % Calculate intersection with fiber's edge
            [ds, medium] = getRefractionPoint(position, direction);
            nIndex = getMediumParams(medium, k);
            
            position(2, :) = position(1, :) + direction*ds;
            
            % Get normal vector of fiber surface
            normalVector = position(2, :);
            normalVector(3) = 0;
            normalVector = normalVector/vecnorm(normalVector);
            
            % Get direction for reflected rays (total or partial)
            reflectedDirection = direction-2*(normalVector*direction')*normalVector;
            
            % Project direction vector on tangent surface and normalize
            tangentVector = direction - (normalVector*direction')*normalVector;
            sineTheta = vecnorm(tangentVector);
            tangentVector = tangentVector/sineTheta;
            
            % Calculate the new angle via Snell's law
            newSineTheta = sineTheta*prevNIndex/nIndex;
            
            if(newSineTheta >= 1)
                % Total internal reflection
                
                % Get reflected direction
                direction = reflectedDirection;
            else
                % Regular refraction
                
                % Get new direction after refraction
                newCosineTheta = sqrt(1-newSineTheta^2);
                newTangentVector = newSineTheta*tangentVector;
                newNormalVector = newCosineTheta*sign(direction*normalVector')*normalVector;
                newDirection = newTangentVector+newNormalVector;

                % Calculate Fresnel coefficients
                cosI = direction*normalVector';
                cosT = newDirection*normalVector';
                fresnelR = (((prevNIndex*cosI-nIndex*cosT)/(prevNIndex*cosI+nIndex*cosT))^2+((prevNIndex*cosT-nIndex*cosI)/(prevNIndex*cosT+nIndex*cosI))^2)/2;
                
                % Random chance to reflect ray
                % Better than splitting into two rays to prevent the power
                % from being diluted
                if rand() < fresnelR
                    newDirection = reflectedDirection;
                end
                
                % Actually change the direction
                direction = newDirection;
            end
            
            % Move in the new direction (out of the boundary)
            position(2, :) = position(2, :) + direction*boundaryDistance;
            ds = ds+boundaryDistance;
            
            medium = getMedium(position(2, :));
        end
        
        % Photon leaves the fiber (and is lost)
        if(medium == 2)
            normalVector = [position(2, 1) ; position(2, 2) ; 0];
            if(direction*normalVector > 0)
                loopOn = false;
                runawayPhotons = runawayPhotons + 1;
            end
        else
            % Photon reaches right end of fiber (and is concentrated)
            if(loopOn && position(2, 3) >= lightL + darkL && position(1, 3) < lightL + darkL)
                loopOn = false;
                
                Pout(k) = Pout(k) + photonPower;
                finalPhotons = finalPhotons + 1;
            end
            
            % Photon reaches ends of fiber but isn't concentrated (and is lost)
            if(loopOn && (position(2, 3) < 0 || position(2, 3) >= lightL + darkL))
                loopOn = false;
                
                runawayPhotons = runawayPhotons + 1;
            end
            
            % A fraction of the power is absorbed by PMMA
            photonPower = photonPower*(1-alfaPMMA*ds);
            
            if photonPower < minimumPower
                loopOn = false;
            end
        end
        
        if medium == 0
            % Index of the N2 interval to use
            j = 1+floor((position(1, 3)+position(2, 3))/2/dz);

            % Photon may be absorbed by dopant (and reemitted)
            if(loopOn && rand() < sigmaabs*(N-N2(j))*ds)
                N2(j) = N2(j)+photonPower*conversionN2(k);
                
                % Generate random direction
                rand1 = 2*pi*rand();
                rand2 = 2*rand()-1;
                direction = [sqrt(1-rand2^2)*cos(rand1) sqrt(1-rand2^2)*sin(rand1) rand2];
                
                % Generate new wavelength after emission
                oldLambda = ll(k);
                k = generateDistributedLambda(ll, emittedDistribution);
                alfaPMMA = alfaPMMAValues(k);
                sigmaabs = sigmaabsValues(k);
                sigmaemi = sigmaemiValues(k);
                
                % Change photon power (Stokes shift loss)
                photonPower = photonPower*oldLambda/ll(k);
                
                % Include non radiative decay chance
                if rand() > quantumYield
                    loopOn = false;
                end
            end
            
            % Photon may cause stimulated emission by the dopant
            if(loopOn && rand() < sigmaemi*N2(j)*ds)
                N2(j) = N2(j) - photonPower*conversionN2(k);
                stimulatedPhotons = stimulatedPhotons + 1;

                runPhoton(position(2, :), direction, k, photonPower);
            end
        end
    end
end

%%%%%%%%%
% Main loop

if nargout == 0
    fprintf('i = 0000000');
end
for i = 1:M
    % Generate random photon from incident sunlight
    k = generateDistributedLambda(ll, solarDistribution);
    photonPower = incomingPower/M;
    
    % Get initial position of photon (directly above fiber, random)
    % and direction (rotated using incidence angle)
    direction = [0 -1 0]*rotationMatrix;
    position = [diameter*(rand()-1/2) 0 lightL*rand()]+[0 (0.5*diameter+10*da)/rotationMatrix(2, 2) 0]*rotationMatrix;

    runPhoton(position, direction, k, photonPower);
    
    if nargout == 0 && mod(i, 5000) == 0
        fprintf('\b\b\b\b\b\b\b%07d', i);
    end
end
if nargout == 0
    fprintf('\n');
end

%%%%%%%%%
% Display results

if nargout == 0
    fprintf('Simulation time: %.1f s\n', toc());
    fprintf('Output power of fiber: %g uW\n', sum(Pout)*1e6);
    fprintf('Photons reaching output: %d/%d\n', finalPhotons, totalPhotons);
    fprintf('Photons leaving the fiber: %d/%d\n', runawayPhotons, totalPhotons);
    fprintf('Stimulated emission photons: %d/%d\n', stimulatedPhotons, totalPhotons);
    
    % Plot power spectrum
    figure(1);
    plot(ll*1e9, Pout*1e-3/dlambda);
    title('Power spectrum at end of fiber (raytracing method)');
    xlabel('\lambda (nm)');
    ylabel('Power spectrum (\muW/nm)');
end

end