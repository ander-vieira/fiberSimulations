function [finalPower] = oneDopantRaytracing(dopant, N, diameter, lightL, darkL, incidenceAngleDegrees)
%ONEDOPANTRAYTRACING Simulate fibers using raytracing as the main tool
%   This function simulates a fiber to obtain a resulting power output by
%   running "photons" through the fiber using raytracing, as opposed to
%   using rate equations

tic;

M = 50000; % Number of photons to simulate

rng('shuffle');

minlambda = 500e-9;
dlambda = 1e-9;
maxlambda = 750e-9;
da = 1e-6;

incidenceAngle = deg2rad(incidenceAngleDegrees);
rotationMatrix = [1 0 0 ; 0 cos(incidenceAngle) -sin(incidenceAngle) ; 0 sin(incidenceAngle) cos(incidenceAngle)]';

[~, sigmaabsFun, sigmaemiFun] = getDyeDopantAttributes(dopant);

ll = minlambda:dlambda:maxlambda;

% Get probability distributions for sunlight photons and emitted photons
solarDistribution = solarIrradianceSpline(ll);
solarConstant = sum(solarDistribution);
solarDistribution = solarDistribution/solarConstant;
solarConstant = solarConstant*dlambda;
emittedDistribution = sigmaemiFun(ll);
emittedDistribution = emittedDistribution/sum(emittedDistribution);

% Initial values
position = zeros(2, 3);
incomingPower = solarConstant*diameter*lightL*cos(incidenceAngle); % W
finalPower = 0; % W
finalPhotons = 0;
runawayPhotons = 0;
PMMAPhotons = 0;
absorptions = 0;

if nargout == 0
    fprintf('i = 000000');
end
for i = 1:M
    % Generate random photon from incident sunlight
    lambda = generateDistributedLambda(ll, solarDistribution);
    alfaPMMA = valuesalfaPMMA(lambda);
    sigmaabs = sigmaabsFun(lambda);
    photonPower = incomingPower/M;
    
    % Get initial position of photon (directly above fiber, random)
    % and direction (rotated using incidence angle)
    direction = [0 -1 0]*rotationMatrix;
    position(2, :) = [diameter*(rand()-1/2) 0 lightL*rand()]+[0 0.55*diameter/rotationMatrix(2, 2) 0]*rotationMatrix;
    
    incoming = true;
    loopOn = true;
    while loopOn
        position(1, :) = position(2, :);
        
        % Get refraction index AND refraction gradient (for GI fibers)
        [prevInFiber, prevN, ~] = getRefractionIndex(position(1, :), diameter, lambda);
        
        % Distance interval travelled 
        ds = prevN*da;
        
        % Change direction based on refraction gradient (for GI fibers
        % only) (see notes)
%         direction = direction + nGradient*da;
%         direction = direction/vecnorm(direction);
        
        % Calculate new position
        position(2, :) = position(1, :) + direction*ds;
        
        [inFiber, n] = getRefractionIndex(position(2, :), diameter, lambda);
        
        if prevInFiber ~= inFiber
            % Refraction in air-fiber interface
            
            incoming = false;
            
            % Get normal vector of fiber surface
            normalVector = position(1, :)+position(2, :);
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
                
                % Keep inside the fiber so photon doesn't disappear
                position(2, :) = position(1, :);
                inFiber = ~inFiber;
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
                
                % Actually change the direction
                direction = newDirection;
            end
        end
        
        % Photon left bounds of simulation (removed for performance)
        if(abs(position(2, 1)) > diameter  || abs(position(2, 2)) > diameter)
            loopOn = false;
        end
        
        % Photon may be absorbed by PMMA (and lost)
        if(inFiber && rand() < alfaPMMA*ds)
            loopOn = false;
            PMMAPhotons = PMMAPhotons + 1;
        end
        
        % Photon may be absorbed by dopant (and reemitted)
        if(inFiber && rand() < sigmaabs*N*ds)
            % Generate random direction
            rand1 = 2*pi*rand();
            rand2 = 2*rand()-1;
            direction = [sqrt(1-rand2^2)*cos(rand1) sqrt(1-rand2^2)*sin(rand1) rand2];
            
            % Generate new wavelength after emission
            oldLambda = lambda;
            lambda = generateDistributedLambda(ll, emittedDistribution);
            alfaPMMA = valuesalfaPMMA(lambda);
            sigmaabs = sigmaabsFun(lambda);
            
            % Change photon power (Stokes shift loss)
            photonPower = photonPower*oldLambda/lambda;
            
            absorptions = absorptions + 1;
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
            finalPower = finalPower + photonPower;
            finalPhotons = finalPhotons + 1;
            loopOn = false;
        end
        
        % Photon reaches ends of fiber but isn't concentrated (and is lost)
        if(~incoming && (position(2, 3) < 0 || position(2, 3) >= lightL + darkL))
            loopOn = false;
        end
    end
    
    if nargout == 0 && mod(i, 100) == 0
        fprintf('\b\b\b\b\b\b%06d', i);
    end
end
if nargout == 0
    fprintf('\n');
end

if nargout == 0
    fprintf('Simulation time: %.1f s\n', toc());
    fprintf('Output power of fiber: %g uW\n', finalPower*1e6);
    fprintf('Photons reaching output: %d/%d\n', finalPhotons, M);
    fprintf('Photons leaving the fiber: %d/%d\n', runawayPhotons, M);
    fprintf('Photons absorbed by PMMA: %d/%d\n', PMMAPhotons, M);
    fprintf('Absorptions: %d\n', absorptions);
end

end

%%%%%%%%%
% Auxiliary functions

% Using a discrete probability distribution, generate random wavelength
function lambda = generateDistributedLambda(ll, distribution)
    randomNumber = rand();
    
    for i = 1:length(ll)
        if(randomNumber < sum(distribution(1:i)))
            break;
        end
    end
    
    lambda = ll(i);
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