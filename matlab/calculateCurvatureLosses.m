function loss = calculateCurvatureLosses(lambda, diameter, curvatureRadius)
%CALCULATECURVATURELOSSES Summary of this function goes here
%   

M = 100000;

ncore = refractionIndexPMMA(lambda);

% Minimum value of cosine of angle so total internal reflection happens
minimumValue = sqrt(1-1/ncore^2);

flatReflectedRays = 0;
reflectedRays = 0;

for i = 1:M
    % Generate random position
    rand1 = 2*pi*rand();
    rand2 = rand()+rand();
    if rand2 > 1
        rand2 = 2 - rand2;
    end
    position = [diameter/2*rand2*cos(rand1) diameter/2*rand2*sin(rand1) 0];

    % Generate random direction
    rand1 = 2*pi*rand();
    rand2 = rand();
    direction = [sqrt(1-rand2^2)*cos(rand1) sqrt(1-rand2^2)*sin(rand1) rand2];

    beta = (direction(1)*position(1)+direction(2)*position(2))/(direction(1)^2+direction(2)^2);
    gamma = (position(1)^2+position(2)^2-diameter^2/4)/(direction(1)^2+direction(2)^2);
    % "Time" until ray reaches edge of fiber
    tauFlat = sqrt(beta^2-gamma)-beta;

    % Get position of first reflection + normal vector in flat fiber
    normalVectorFlat = [position(1)+direction(1)*tauFlat position(2)+direction(2)*tauFlat 0];
    normalVectorFlat = normalVectorFlat/vecnorm(normalVectorFlat);
    
    % Reflection angle
    cosAngleFlat = normalVectorFlat*direction';
    
    if cosAngleFlat < minimumValue
        flatReflectedRays = flatReflectedRays+1;
    end
    
    vectorP = position - [curvatureRadius 0 0];
    vectorV = direction';
    vectorPflat = [vectorP(1) 0 vectorP(3)];
    vectorVflat = [vectorV(1) 0 vectorV(3)]';
    normP = vecnorm(vectorP);
    normV = vecnorm(vectorV);
    normPflat = vecnorm(vectorPflat);
    normVflat = vecnorm(vectorVflat);
    PV = vectorP*vectorV;
    PVflat = vectorPflat*vectorVflat;
    
    polynom = zeros(1, 5);
    polynom(1) = normV^4;
    polynom(2) = 4*normV^2*PV;
    polynom(3) = 4*PV^2+2*normV^2*(normP^2+curvatureRadius^2-diameter^2/4)-4*curvatureRadius^2*normVflat^2;
    polynom(4) = 4*PV*(normP^2+curvatureRadius^2-diameter^2/4)-8*curvatureRadius^2*PVflat;
    polynom(5) = (normP^2+curvatureRadius^2-diameter^2/4)^2-4*curvatureRadius^2*normPflat^2;
    
    tauValues = roots(polynom);
    tau = min(tauValues(tauValues==real(tauValues) & tauValues>0));
    
    newPosition = position+tau*direction;
    
    forwardVector = [newPosition(3) 0 curvatureRadius-newPosition(1)];
    forwardVector = forwardVector/vecnorm(forwardVector);
    
    normalVector = newPosition-(newPosition*forwardVector')*forwardVector;
    normalVector = normalVector/vecnorm(normalVector);
    
    % Reflection angle
    cosAngle = normalVector*direction';
    
    % If total internal reflection here, it will always get TIR because of
    % symmetry
    if cosAngle < minimumValue
        reflectedRays = reflectedRays+1;
    end
end

loss = reflectedRays/flatReflectedRays;

end

