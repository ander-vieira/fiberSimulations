function [finalPower] = oneDopantRaytracing(dopant, N, diameter, lightL, darkL, incidenceAngleDegrees)
%ONEDOPANTRAYTRACING Simulate fibers using raytracing as the main tool
%   This function simulates a fiber to obtain a resulting power output by
%   running "photons" through the fiber using raytracing, as opposed to
%   using rate equations

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

minlambda = 250e-9;
dlambda = 2e-9;
maxlambda = 750e-9;
ll = minlambda:dlambda:maxlambda;

da = 5e-5;

% Rotation due to sunlight incidence angle
incidenceAngle = deg2rad(incidenceAngleDegrees);
rotationMatrix = [1 0 0 ; 0 cos(incidenceAngle) -sin(incidenceAngle) ; 0 sin(incidenceAngle) cos(incidenceAngle)]';

[tauRad, sigmaabsFun, sigmaemiFun, tauNR] = getDyeDopantAttributes(dopant);

% Get probability distributions for sunlight photons and emitted photons
solarDistribution = solarIrradianceSpline(ll);
solarConstant = sum(solarDistribution);
solarDistribution = solarDistribution/solarConstant;
solarConstant = solarConstant*dlambda; % W/m^2
sigmaabsValues = sigmaabsFun(ll);
sigmaemiValues = sigmaemiFun(ll);
emittedDistribution = sigmaemiValues/sum(sigmaemiValues);
alfaPMMAValues = attenuationPMMA(ll);

quantumYield = tauNR/(tauRad+tauNR);
% Conversion constant from power to N2
conversionN2 = 4*ll/(pi*h*c*diameter^2*dz)/(1/tauRad+1/tauNR); % m^-3/W

% Initial values
N2 = zeros(1, numzz);
incomingPower = solarConstant*diameter*lightL*cos(incidenceAngle); % W
minimumPower = incomingPower/M/100;
finalPower = 0; % W
totalPhotons = 0;
finalPhotons = 0;
runawayPhotons = 0;
PMMAPhotons = 0;
stimulatedPhotons = 0;

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

        % Get refraction index AND refraction gradient (for GI fibers)
        [prevInFiber, prevN, ~] = getRefractionIndex(position(1, :), diameter, ll(k));
        
        % Distance interval travelled 
        ds = prevN*da;
        
        % Change direction based on refraction gradient (for GI fibers
        % only) (see notes)
%         direction = direction + nGradient*da;
%         direction = direction/vecnorm(direction);

        % Calculate new position
        position(2, :) = position(1, :) + direction*ds;
        
        [inFiber, n] = getRefractionIndex(position(2, :), diameter, ll(k));
        
        if prevInFiber ~= inFiber
            % Refraction in air-fiber interface
            
            % Calculate intersection with fiber's edge
            beta = (direction(1)*position(1, 1)+direction(2)*position(1, 2))/(direction(1)^2+direction(2)^2);
            gamma = (position(1, 1)^2+position(1, 2)^2-diameter^2/4)/(direction(1)^2+direction(2)^2);
            
            if beta > 0
                ds = sqrt(beta^2-gamma)-beta;
            else
                ds = -sqrt(beta^2-gamma)-beta;
            end
            position(2, :) = position(1, :) + direction*ds;
            
            % Get normal vector of fiber surface
            normalVector = position(2, :);
            normalVector(3) = 0;
            normalVector = normalVector/vecnorm(normalVector);
            
            % Projection of direction vector on fiber surface
            projectedVector = direction - (normalVector*direction')*normalVector;
            
            % Scale projected vector via Snell's law
            newProjectedVector = projectedVector*prevN/n;
            
            if(vecnorm(newProjectedVector) >= 1)
                % Total internal reflection
                
                % Get reflected direction
                direction = 2*projectedVector - direction;
            else
                % Regular refraction
                
                % Get new direction after refraction
                newDirection = newProjectedVector + sign(direction*normalVector')*sqrt(1 - vecnorm(newProjectedVector)^2)*normalVector;

                % Calculate Fresnel coefficients
                cosI = direction*normalVector';
                cosT = newDirection*normalVector';
                fresnelR = (((prevN*cosI-n*cosT)/(prevN*cosI+n*cosT))^2+((prevN*cosT-n*cosI)/(prevN*cosT+n*cosI))^2)/2;
                fresnelT = 1-fresnelR;
                
                % Refracted power kept in current photon
                photonPower = photonPower*fresnelT;
                
                % TODO: generate new photon with reflected part of power
                reflectedDirection = 2*projectedVector - direction;
                reflectedPower = photonPower*fresnelR;
                if reflectedPower > minimumPower
                    runPhoton(position(2, :)+reflectedDirection*prevN*da, reflectedDirection, k, reflectedPower);
                end
                
                % Actually change the direction
                direction = newDirection;
            end
            
            % Move in the new direction (out of the boundary)
            position(2, :) = position(2, :) + direction*prevN*da;
            ds = ds + prevN*da;
            
            [inFiber, ~] = getRefractionIndex(position(2, :), diameter, ll(k));
        end
        
        % Photon leaves the fiber (and is lost)
        if(~inFiber)
            normalVector = [position(2, 1) ; position(2, 2) ; 0];
            if(direction*normalVector > 0)
                loopOn = false;
                runawayPhotons = runawayPhotons + 1;
            end
        end
        
        % Photon reaches right end of fiber (and is concentrated)
        if(loopOn && inFiber && position(2, 3) >= lightL + darkL && position(1, 3) < lightL + darkL)
            loopOn = false;
            
            finalPower = finalPower + photonPower;
            finalPhotons = finalPhotons + 1;
        end
        
        % Photon reaches ends of fiber but isn't concentrated (and is lost)
        if(loopOn && inFiber && (position(2, 3) < 0 || position(2, 3) >= lightL + darkL))
            loopOn = false;
            
            runawayPhotons = runawayPhotons + 1;
        end
        
        % Photon may be absorbed by PMMA (and lost)
        if(inFiber && rand() < alfaPMMA*ds)
            loopOn = false;
            
            PMMAPhotons = PMMAPhotons + 1;
        end
        
        % Index of the N2 interval to use
        j = 1+floor((position(1, 3)+position(2, 3))/2/dz);
        
        % Photon may be absorbed by dopant (and reemitted)
        if(loopOn && inFiber && rand() < sigmaabs*(N-N2(j))*ds)
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
        if(loopOn && inFiber && rand() < sigmaemi*N2(j)*ds)
            N2(j) = N2(j) - photonPower*conversionN2(k);
            stimulatedPhotons = stimulatedPhotons + 1;
            
            runPhoton(position(2, :), direction, k, photonPower);
        end
    end
end

if nargout == 0
    fprintf('i = 000000');
end
for i = 1:M
    % Generate random photon from incident sunlight
    k = generateDistributedLambda(ll, solarDistribution);
%     alfaPMMA = alfaPMMAValues(k);
%     sigmaabs = sigmaabsValues(k);
    photonPower = incomingPower/M;
    
    % Get initial position of photon (directly above fiber, random)
    % and direction (rotated using incidence angle)
    direction = [0 -1 0]*rotationMatrix;
%     position(2, :) = [diameter*(rand()-1/2) 0 lightL*rand()]+[0 0.55*diameter/rotationMatrix(2, 2) 0]*rotationMatrix;
      position = [diameter*(rand()-1/2) 0 lightL*rand()]+[0 0.55*diameter/rotationMatrix(2, 2) 0]*rotationMatrix;

    runPhoton(position, direction, k, photonPower);
    
    if nargout == 0 && mod(i, 5000) == 0
        fprintf('\b\b\b\b\b\b%06d', i);
    end
end
if nargout == 0
    fprintf('\n');
end

if nargout == 0
    fprintf('Simulation time: %.1f s\n', toc());
    fprintf('Output power of fiber: %g uW\n', finalPower*1e6);
    fprintf('Photons reaching output: %d/%d\n', finalPhotons, totalPhotons);
    fprintf('Photons leaving the fiber: %d/%d\n', runawayPhotons, totalPhotons);
    fprintf('Photons absorbed by PMMA: %d/%d\n', PMMAPhotons, totalPhotons);
    fprintf('Stimulated emission photons: %d/%d\n', stimulatedPhotons, totalPhotons);
end

end

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

% Get refraction index and refraction gradient for a position
% Needs to be defined based on fiber/system geometry
% Gradient only matters for GI fibers
function [inFiber, n, nGradient] = getRefractionIndex(position, diameter, lambda)
    distanceToCenter = position(1)^2+position(2)^2;
    
    if(distanceToCenter < diameter^2/4)
        % Position inside fiber
    
        inFiber = true;
        n = refractionIndexPMMA(lambda);
        % This is for gradient-index fibers, ignore for now
        nGradient = 0;
    else
        % Position outside fiber
        
        inFiber = false;
        n = 1;
        % This is for gradient-index fibers, ignore for now
        nGradient = 0;
    end
end