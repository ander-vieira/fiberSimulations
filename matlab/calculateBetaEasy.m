function beta = calculateBetaEasy(lambda)
%CALCULATEBETAEASY Summary of this function goes here
%   Detailed explanation goes here

M = 1000000;

ncore = refractionIndexPMMA(lambda);

% Minimum value of cosine of angle so total internal reflection happens
minimumValue = sqrt(1-1/ncore^2);

reflectedRays = 0;

for i=1:M
    % Generate random position
    rand1 = 2*pi*rand();
    rand2 = rand()+rand();
    if rand2 > 1
        rand2 = 2 - rand2;
    end
    position = [rand2*cos(rand1) rand2*sin(rand1) 0];
    
    % Generate random direction
    rand1 = 2*pi*rand();
    rand2 = rand();
    direction = [sqrt(1-rand2^2)*cos(rand1) sqrt(1-rand2^2)*sin(rand1) rand2];
    
    beta = (direction(1)*position(1)+direction(2)*position(2))/(direction(1)^2+direction(2)^2);
    gamma = (position(1)^2+position(2)^2-1)/(direction(1)^2+direction(2)^2);
    
    % This can't happen if the ray comes from inside the fiber
    if gamma>1 || beta^2-gamma < 0
        fprintf('How?\n');
    end
    
    % "Time" until ray reaches edge of fiber
    tau = sqrt(beta^2-gamma)-beta;
    
    % Get position of first reflection + normal vector
    normalVector = [position(1)+direction(1)*tau position(2)+direction(2)*tau 0];
    
    % Reflection angle
    cosAngle = normalVector*direction';
    
    % If total internal reflection here, it will always get TIR because of
    % symmetry
    if cosAngle < minimumValue
        reflectedRays = reflectedRays+1;
    end
end

beta = reflectedRays/M/2;

end

