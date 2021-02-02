function [wnsp] = oneDopantReverseW(dopant, N, diameter, lightL, darkL, Poutfun)
%ONEDOPANTREVERSEW Reverse red-shift in a fiber's emission spectrum
%   By applying the reverse process of that normally used in the
%   iterative simulation method, you can extract the shape of the
%   fluorescence spectrum of the dopant from the red-shifted one.
%
%   Parameters:
%   dopant: name of the dye dopant to process
%   N: concentration of the dopant (m^-3)
%   diameter: diameter of the fiber (m)
%   lightL: length of the fiber illuminated by sunlight
%   darkL: length of the fiber not illuminated at the end
%   Poutfun: a function that returns the output power as measured, as a
%   function of lambda (wavelength)

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

Pout = Poutfun(ll);

alfaPMMA = valuesalfaPMMA(ll);

isol = solarIrradianceSpline(ll);

ncore = refractionIndexPMMA(ll);
% beta = (ncore - 1)./(2*ncore);
beta = zeros(1, numll);

efficiency = zeros(1, numll);
for k = 1:numll
    alfaCore = alfaPMMA(k) + sigmaabs(k)*N + realmin;
%     efficiency(k) = fiberAbsorptionNoReflections(ncore(k), diameter, sigmaabs(k)*N, alfaCore);
%     efficiency(k) = fiberAbsorptionReflections(ncore(k), diameter, sigmaabs(k)*N, alfaCore);
    efficiency(k) = fiberAbsorptionTwoInterfaces(ncore(k), 1.4, diameter, .98, sigmaabs(k)*N, alfaCore, alfaPMMA(k));
    beta(k) = calculateBetaIntegral(ll(k));
end

% Precalculated constants
concentrationToPower = pi*h*c*diameter^2./(4*ll);
Nsolconst = diameter*sum(isol*dlambda.*efficiency./concentrationToPower);
Nabsconst = sigmaabs./concentrationToPower;
Pattconst = (alfaPMMA+N*sigmaabs)*dz;
PNconst1 = concentrationToPower.*beta.*wnsp*dz/tau;

% Outer loop: in each iteration, refine the values of wnsp so they match
% the output power better
errorOuter = 1;
while errorOuter > 1e-5
    oldWnsp = wnsp;
    
    % Scale up the new spectrum using F-L equation
    sigmaemiScale = fuchtbauerLadenburgDiscrete(tau, dlambda, ll, wnsp);
    sigmaemi = wnsp*sigmaemiScale;
    
    % These constants need to be recalculated since they depend on the
    % emission spectrum
    Nestconst = sigmaemi./concentrationToPower;
    PNconst2 = (sigmaabs+sigmaemi)*dz;

    % This part is the same as in oneDopantIterative: regular simulation
    % of the fiber
    P = zeros(numzz, numll);
    Pleft = zeros(numzz, numll);
    N2 = zeros(numzz, 1);

    errorInner = 1;
    while errorInner > 1e-8
        % Boundary condition for P
        P(1, :) = zeros(1, numll);
        Pleft(end, :) = zeros(1, numll);

        previousP = P(end, :);

        % Update N2
        for j = 1:numzz-1
            evalP = (P(j, :)+P(j+1, :)+Pleft(j, :)+Pleft(j+1, :))/2;
            wabs = sum(Nabsconst.*evalP);
            west = sum(Nestconst.*evalP);

            A = 1/tau+wabs+west;

            if j <= lightj
                b = Nsolconst+N*wabs;
            else
                b = N*wabs;
            end

            N2(j) = A\b;
        end

        % Update P
        for j = 1:numzz-1
            evalN2 = N2(j);

            for k = 1:numll
                evalP = P(j, k);

                P(j+1, k) = P(j, k);

                P(j+1, k) = P(j+1, k) - Pattconst(k)*evalP;
                P(j+1, k) = P(j+1, k) + PNconst1(k)*evalN2;
                P(j+1, k) = P(j+1, k) + PNconst2(k)*evalN2*evalP;
            end
        end

        % Update Pleft
        for jinv = numzz:-1:2
            evalN2 = N2(jinv-1);

            for k = 1:numll
                evalP = Pleft(jinv, k);

                Pleft(jinv-1, k) = Pleft(jinv, k);

                Pleft(jinv-1, k) = Pleft(jinv-1, k) - Pattconst(k)*evalP;
                Pleft(jinv-1, k) = Pleft(jinv-1, k) + PNconst1(k)*evalN2;
                Pleft(jinv-1, k) = Pleft(jinv-1, k) + PNconst2(k)*evalN2*evalP;
            end
        end

        % Calculate the error (for the loop condition)
        errorInner = max(abs((P(end, :)-previousP)./(P(end, :)+realmin)));
    end
    
    % After getting P and N2, you also need the derivatives of P over each
    % individual weight (each component of the wnsp vector)
    % To calculate it you differentiate the rate equations and get new
    % rate equations for the derivatives (called Q here)
    % You only need values of Q for the end of the fiber
    % See notes for full derivation
    Qend = zeros(numll, numll);
    
    % For each wavelength (each component of wnsp) you need to iterate
    % The method is similar to the one in oneDopantIterative
    % But the equations to solve are different
    for q = 1:numll
        % Initialize
        Q = zeros(numzz, numll);
        Tn = zeros(1, numzz);
        
        errorInner = 1;
        while errorInner > 1e-8
            % Boundary condition for Q
            % Only propagate Q in one direction: you could do both in
            % theory, but the derivatives don't need that much precision
            Q(1, :) = zeros(1, numll);

            previousQ = Q(end, :);

            % Update TN
            for j = 1:numzz-1
                evalP = (P(j, :)+P(j+1, :)+Pleft(j, :)+Pleft(j+1, :))/2;
                evalQ = (Q(j, :)+Q(j+1, :))/2;
                evalN2 = N2(j);
                wabs = sum(Nabsconst.*evalP);
                west = sum(Nestconst.*evalP);
                Yabs = sum(Nabsconst.*evalQ);
                Yest = sum(Nestconst.*evalQ)+(sigmaemi(q)+sigmaemiScale)*evalP(q)/concentrationToPower(q)-west;

                A = 1/tau+wabs+west;
                
                b = Yabs*N - (Yabs + Yest)*evalN2;

                Tn(j) = A\b;
            end

            % Update Q
            for j = 1:numzz-1
                evalN2 = N2(j);
                evalTn = Tn(j);

                for k = 1:numll
                    evalP = P(j, k);
                    evalQ = Q(j, k);

                    Q(j+1, k) = Q(j, k);

                    % These terms are equivalent to the ones from
                    % oneDopantIterative
                    Q(j+1, k) = Q(j+1, k) - Pattconst(k)*evalQ;
                    Q(j+1, k) = Q(j+1, k) + PNconst1(k)*evalTn;
                    Q(j+1, k) = Q(j+1, k) + PNconst2(k)*evalN2*evalQ;
                    
                    % These new terms appear when differentiating
                    Q(j+1, k) = Q(j+1, k) + PNconst2(k)*evalTn*evalP;
                    
                    % The DIAGONAL terms in Q have different behaviour
                    if k == q
                        Q(j+1, k) = Q(j+1, k) + sigmaemiScale*dz*evalN2*evalP;
                        Q(j+1, k) = Q(j+1, k) + concentrationToPower(k)*beta(k)*dz*evalN2;
                    else
                        Q(j+1, k) = Q(j+1, k) - sigmaemi(k)*dz*evalN2*evalP;
                        Q(j+1, k) = Q(j+1, k) - PNconst1(k)*evalN2;
                    end
                end
            end

            % Calculate the error (for the loop condition)
            errorInner = max(abs((Q(end, :)-previousQ)./(Q(end, :)+realmin)));
        end
        
        Qend(:, q) = Q(end, :)';
    end
    
    % Get the difference between expected and simulated powers
    dif = Pout - P(end, :);
    
    % Applying N-R in matrix form, get new values for the wnsp
    wnsp = wnsp - Qend\dif';

    errorOuter = max(abs(wnsp-oldWnsp));
end

end

% Apply the F-L equation to a discrete sample of values
% Returns the value the wnsp needs to be scaled by
function scale = fuchtbauerLadenburgDiscrete(tau, dlambda, ll, wnsp)

% Constants
c = 3e8;% m/s

ncore = refractionIndexPMMA(ll);

% Replace the integral with a sum since it's discrete
integralValue = sum(wnsp.*ncore.^2*dlambda./ll.^4);

scale = 1/(8*pi*c*tau*integralValue);

end