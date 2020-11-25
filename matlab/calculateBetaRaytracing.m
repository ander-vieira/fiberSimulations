function beta = calculateBetaRaytracing(diameter, L, lambda)
%CALCULATEBETARAYTRACING Summary of this function goes here
%   Detailed explanation goes here

tic;

M = 10000; % Number of photons to simulate

rng('shuffle');

da = 1e-5;

position = zeros(2, 3);
rightPhotons = 0;
lostLength = zeros(1, M);

if nargout == 0
    fprintf('i = 000000');
end
for i = 1:M
    % Generate random position
    rand1 = 2*pi*rand();
    rand2 = rand()+rand();
    if rand2 > 1
        rand2 = 2 - rand2;
    end
    position(2, :) = [diameter/2*rand2*cos(rand1) diameter/2*rand2*sin(rand1) 0];
    
    % Generate random direction
    rand2 = 0;
    while rand2 < 0.001 % These would just stay in the fiber forever
        rand1 = 2*pi*rand();
        rand2 = rand();
        direction = [sqrt(1-rand2^2)*cos(rand1) sqrt(1-rand2^2)*sin(rand1) rand2];
    end
    
    loopOn = true;
    while loopOn
        position(1, :) = position(2, :);
        
        % Get refraction index AND refraction gradient (for GI fibers)
        [~, prevN] = getRefractionIndex(position(1, :), diameter, lambda);
        
        % Distance interval travelled 
        ds = prevN*da;
        
        % Change direction based on refraction gradient (for GI fibers
        % only) (see notes)
%         direction = direction + nGradient*da;
%         direction = direction/vecnorm(direction);

        % Calculate new position
        position(2, :) = position(1, :) + direction*ds;
        
        [inFiber, n] = getRefractionIndex(position(2, :), diameter, lambda);
        
        if ~inFiber
            % Refraction in air-fiber interface
            
            % Calculate intersection with fiber's edge
            beta = (direction(1)*position(1, 1)+direction(2)*position(1, 2))/(direction(1)^2+direction(2)^2);
            gamma = (position(1, 1)^2+position(1, 2)^2-diameter^2/4)/(direction(1)^2+direction(2)^2);
            
            ds = sqrt(beta^2-gamma)-beta;
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
                
                % Keep inside the fiber so photon doesn't disappear
%                 position(2, :) = position(1, :);
            else
                % Regular refraction
                lostLength(i) = position(2, 3);
                
                loopOn = false;
            end
        end
        
        % Photon reaches left end of fiber (and is counted)
        if(loopOn && position(2, 3) < 0)
            loopOn = false;
%             leftPhotons = leftPhotons + 1;
        end
        
        % Photon reaches right end of fiber (and is counted)
        if(loopOn && position(2, 3) >= L)
            loopOn = false;
            rightPhotons = rightPhotons + 1;
            lostLength(i) = L;
        end
    end
    if nargout == 0 && mod(i, 100) == 0
        fprintf('\b\b\b\b\b\b%06d', i);
    end
end
if nargout == 0
    fprintf('\n');
end

beta = rightPhotons/2/M;

if nargout == 0
    fprintf('Simulation time: %.1f s\n', toc());
    fprintf('Average beta: %g\n', beta);
    
    figure(1);
    histogram(lostLength*1e2, 200);
    xlabel('Length (cm)');
end

end

%%%%%%%%%
% Auxiliary functions

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
