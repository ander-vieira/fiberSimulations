function [beta, KZCore, KZClad] = geometricalParams2(ncore, nclad, q)
% GEOMETRICALPARAMS2
%

% Limits for the cosines of the angles of incidence into each interphase
% If cos < M, then sin > sinCut, and TIR happens in that interphase
M1 = sqrt(1-(nclad/ncore)^2);
M2 = sqrt(1-(1/nclad)^2);

% For a given ray defined by u0 and alpha (see notes), calculate parameters
% TIR: whether TIR happens to keep the ray in the fiber (1 or 0)
% F: fraction of distance traversed inside the core (between 0 and 1)
function [TIR, invKZCore, invKZClad] = solveRayParams(u0, alpha)
    % Initial position and velocity of the ray
    p0 = [q*u0 0 0];
    v1 = [0 sin(alpha) cos(alpha)];
    
    % Intersection of the ray with core-cladding interphase
    p1 = [q*u0 q*sqrt(1-u0.^2) q*sqrt(1-u0.^2)/tan(alpha)];
    
    % Distance traveled in the core
    distance1 = p1-p0;
    
    % Check TIR in the core-cladding interphase
    % Normal vector is already normalized by definition
    normalVector1 = p1;
    normalVector1(3) = 0;
    cosTheta1 = v1*normalVector1';
    
    if cosTheta1 < M1
        % If TIR in the core-cladding interphase: stays inside the core
        TIR = 1;
        
        invKZCore = distance1(3)/vecnorm(distance1);
        invKZClad = 0;
    else
        % If refraction: goes on to cladding, possibly TIR in the next
        % interphase
        
        % Calculate velocity vector after refraction via Snell's law
        normalV1 = cosTheta1*normalVector1;
        tangentV1 = v1-normalV1;
        tangentV2 = tangentV1*ncore/nclad;
        normalV2 = sqrt(cosTheta1^2-1+(nclad/ncore)^2)*normalVector1;
        v2 = tangentV2+normalV2;
        
        % Calculate intersection with cladding-air interphase
        A = (p1(1)*v2(1)+p1(2)*v2(2))/(v2(1)^2+v2(2)^2);
        B = (1-p1(1)^2-p1(2)^2)/(v2(1)^2+v2(2)^2);
        tau = sqrt(A^2+B)-A;
        
        % Intersection point with cladding-air interphase
        p2 = p1+v2*tau;
        
        % Check TIR in the core-cladding interphase
        % Normal vector is already normalized by definition
        normalVector2 = p2;
        normalVector2(3) = 0;
        cosTheta2 = v2*normalVector2';
        
        % Distance traveled in the cladding
        distance2 = p2-p1;
        
        if cosTheta2 < M2
            % If TIR in the cladding-air interphase
            TIR = 1;
            
            % Calculate the inverses of KZcore and KZclad, which measure
            % the fraction between distance along the z axis and distance
            % traversed (including zigzag movement) in core and cladding,
            % respectively
            invKZCore = (distance1(3)+distance2(3))/vecnorm(distance1);
            invKZClad = (distance1(3)+distance2(3))/vecnorm(distance2);
        else
            % If no TIR, the ray is lost, other values aren't used
            TIR = 0;
            
            invKZCore = 0;
            invKZClad = 0;
        end
    end
end

% Coding G-K quadrature by hand is far more efficient in this case
% due to the multiple integrals that can be computed at the same time,
% plus evaluations of solveRayParams are very costly.

% Coefficients for G-K quadrature (order 15)
gaussKronrodX = [-0.99145 -0.94911 -0.86486 -0.74153 -0.58609 -0.40585 -0.20778 0 0.20778 0.40585 0.58609 0.74153 0.86486 0.94911 0.99145];
gaussKronrodW = [0.02294 0.06309 0.10479 0.14065 0.16900 0.19035 0.20443 0.20948 0.20443 0.19035 0.16900 0.14065 0.10479 0.06309 0.02294];

% Coefficients for G-K quadrature (order 29)
% gaussKronrodX = [-0.9977205937565431220162
% -0.9862838086968123388416
% -0.9631583382788532040077
% -0.9284348836635735173364
% -0.8829146632520570399352
% -0.8272013150697649931898
% -0.7617567525622055250577
% -0.687292904811685470148
% -0.604789365940921592412
% -0.5152486363581540919653
% -0.4196558976429790001574
% -0.3191123689278897604357
% -0.2148359185334849021051
% -0.1080549487073436620662
% 0
% 0.1080549487073436620662
% 0.214835918533484902105
% 0.3191123689278897604357
% 0.4196558976429790001574
% 0.5152486363581540919653
% 0.604789365940921592412
% 0.687292904811685470148
% 0.7617567525622055250577
% 0.8272013150697649931898
% 0.8829146632520570399352
% 0.9284348836635735173364
% 0.9631583382788532040077
% 0.9862838086968123388416
% 0.9977205937565431220162];
% gaussKronrodW = [0.00613955868637813137607
% 0.017148458909935506453
% 0.029048701261508505535
% 0.040250594872688612167
% 0.05069154326046537607
% 0.06066712586674214967355
% 0.0701029790027469826632
% 0.07865579724962168225
% 0.0861837628694581380093
% 0.092736830017853571435
% 0.098264926472103737472
% 0.1026166273213997903287
% 0.105731639841764335702
% 0.1076264211141186523633
% 0.1082700665064296570044
% 0.1076264211141186523633
% 0.1057316398417643357017
% 0.1026166273213997903287
% 0.098264926472103737472
% 0.0927368300178535714355
% 0.0861837628694581380093
% 0.07865579724962168225
% 0.070102979002746982663
% 0.0606671258667421496735
% 0.05069154326046537607
% 0.0402505948726886121668
% 0.029048701261508505535
% 0.017148458909935506453
% 0.0061395586863781313761];

% Normalize G-K weights to 1 (should normalize to 2 but it cancels out)
gaussKronrodW = gaussKronrodW/sum(gaussKronrodW);

% Integrals are calculated by adding components
integral1 = 0;
integral2 = 0;
integral3 = 0;
integral4 = 0;
% Iterate over all G-K nodes for each variable
for i = 1:length(gaussKronrodX)
    u0 = (1+gaussKronrodX(i))/2;
    for j = 1:length(gaussKronrodX)
        alpha = (1+gaussKronrodX(j))*pi/4;
        
        % Evaluate solveRayParams
        [TIR, invKZCore, invKZClad] = solveRayParams(u0, alpha);
        
        % Evaluate the integrands use solveRayParams results
        % and add each integrand evaluation to the respective integral sum
        integrand1 = sqrt(1-u0^2)*sin(alpha);
        integral1 = integral1 + pi/2*gaussKronrodW(i)*gaussKronrodW(j)*integrand1;
        integrand2 = integrand1*TIR;
        integral2 = integral2 + pi/2*gaussKronrodW(i)*gaussKronrodW(j)*integrand2;
        integrand3 = integrand2*invKZCore;
        integral3 = integral3 + pi/2*gaussKronrodW(i)*gaussKronrodW(j)*integrand3;
        integrand4 = integrand2*invKZClad;
        integral4 = integral4 + pi/2*gaussKronrodW(i)*gaussKronrodW(j)*integrand4;
    end
end

% Using the computed integrals, get the values of the parameters
beta = integral2/integral1/2;
KZCore = integral2/integral3;
KZClad = integral2/integral4;

% % Slow method using built-in integral2 method
% % Gives very similar results, depending on order of G-K used
% integrand1 = @(u0, alpha) sin(alpha).*sqrt(1-u0.^2);
% 
% function result = integrand2(u0, alpha)
%     if length(u0) > 1
%         result = zeros(size(u0));
%         for i=1:size(u0, 1)
%             for j=1:size(u0, 2)
%                 result(i, j) = integrand2(u0(i, j), alpha(i, j));
%             end
%         end
%     else
%         [TIR, ~] = solveRayParams(u0, alpha);
% 
%         result = integrand1(u0, alpha)*TIR;
%     end
% end
% 
% function result = integrand3(u0, alpha)
%     if length(u0) > 1
%         result = zeros(size(u0));
%         for i=1:size(u0, 1)
%             for j=1:size(u0, 2)
%                 result(i, j) = integrand3(u0(i, j), alpha(i, j));
%             end
%         end
%     else
%         [TIR, invKZCore] = solveRayParams(u0, alpha);
% 
%         result = integrand1(u0, alpha)*TIR*invKZCore;
%     end
% end
% 
% function result = integrand4(u0, alpha)
%     if length(u0) > 1
%         result = zeros(size(u0));
%         for i=1:size(u0, 1)
%             for j=1:size(u0, 2)
%                 result(i, j) = integrand4(u0(i, j), alpha(i, j));
%             end
%         end
%     else
%         [TIR, ~, invKZClad] = solveRayParams(u0, alpha);
% 
%         result = integrand1(u0, alpha)*TIR*invKZClad;
%     end
% end
% 
% result1 = integral2(integrand1, 0, 1, 0, pi/2);
% result2 = integral2(@integrand2, 0, 1, 0, pi/2);
% result3 = integral2(@integrand3, 0, 1, 0, pi/2);
% result4 = integral2(@integrand4, 0, 1, 0, pi/2);
% beta = result2/result1/2;
% KZCore = result2/result3;
% KZClad = result2/result4;

end