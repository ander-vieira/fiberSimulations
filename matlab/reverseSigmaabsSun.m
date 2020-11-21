function [sigmaabs] = reverseSigmaabsSun(N, diameter, lightL, darkL, ll, wnsp, Pout)
%REVERSESIGMAABSSUN Summary of this function goes here
%   Detailed explanation goes here

tic;

ncore = refractionIndexPMMA(ll);
alfaPMMA = valuesalfaPMMA(ll);
isol = solarIrradianceSpline(ll, "surface");

beta = (ncore - 1)./(2*ncore);

% Interval used for numerical differentiation
dsigma = 1e-24;

% Iterate using Newton-Raphson adapted for M-sized vectors
sigmaabs = zeros(1, length(ll));
error = 1;
while error > 1e-5
    % Get expected power spectrum shape
    [fun, J] = solveSystem(ll, N, diameter, lightL, darkL, sigmaabs, alfaPMMA, isol, beta, wnsp, dsigma);
    
    deltaS = (J\(Pout'-fun))';
    
    % Calculate relative error in each step
    error = max(abs(deltaS./sigmaabs));
    
    % Calculate new values for the absorption sigma
    sigmaabs = sigmaabs + deltaS;
    
    sigmaabs = abs(sigmaabs)+realmin;
end
    
fprintf("Time elapsed in reverseSigmaabsSun: %.2f\n", toc());

end

% Simulates the behaviour of the POF via the stationary method
% Returns the shape of the output power spectrum
% Also returns the numerical differential of the shape of P
% with respect to each value in WNSP for N-R purposes
function [result, diff] = solveSystem(ll, N, diameter, lightL, darkL, sigmaabs, alfaPMMA, isol, beta, wnsp, dsigma)
    M = length(ll);
    diff = zeros(M);

%     efficiency = zeros(1, M);
    efficiency = (N*sigmaabs)./(N*sigmaabs+alfaPMMA+realmin).*(1-exp(-alfaPMMA-N*sigmaabs*diameter*0.919));
%     for k = 1:M
%         alfaCore = alfaPMMA(k) + N*sigmaabs(k) + realmin;
%         efficiency(k) = fiberAbsorptionTwoInterfaces(ncore(k), 1.4, diameter, .98, sigmaabs(k)*N, alfaCore, alfaPMMA(k));
%     end
    
    % Define C and G, which constitute a system of linear differential
    % equations
    C = sum(isol.*efficiency.*ll*2e-9)*diameter*beta.*wnsp./ll;

    G = zeros(M);
    Gprime = zeros(M);
    for i = 1:M
        for j = 1:M
            if(i == j)
                G(i, i) = -sigmaabs(i).*N.*(1-beta(i).*wnsp(i))-alfaPMMA(i)+sqrt(realmin);
            else
                G(j, i) = sigmaabs(i)*N*beta(j)*wnsp(j)*ll(i)/ll(j);
            end
        end
    end
    
    % Get the shape of the power output
    % via solving a linear differential equation system
    % defined by G and C
    result = expm(G*darkL)*(expm(G*lightL)-eye(M))*(G\C');
    
    % Get the differential for each component k of WNSP
    for k = 1:M
        % Increase component k of WNSP slightly
        direction = [zeros(1,k-1) 1 zeros(1,M-k)];
        newSigma = (sigmaabs+dsigma*direction);
        
        efficiency = (N*newSigma)./(N*newSigma+alfaPMMA+realmin).*(1-exp(-alfaPMMA-N*newSigma*diameter*0.919));
    
        % Define a new C and G for the modified WNSP
        Cprime = sum(isol.*efficiency.*ll*2e-9)*diameter*beta.*wnsp./ll;
    
        for i = 1:M
            for j = 1:M
                if(i == j)
                    Gprime(i, i) = -newSigma(i).*N.*(1-beta(i).*wnsp(i))-alfaPMMA(i)+sqrt(realmin);
                else
                    Gprime(j, i) = newSigma(i)*N*beta(j)*wnsp(j)*ll(i)/ll(j);
                end
            end
        end
        
        % Solve the system of linear differential equations again
        resultprime = expm(Gprime*darkL)*(expm(Gprime*lightL)-eye(M))*(Gprime\Cprime');
        
        % Get the derivative from substracting both results
        diff(:, k) = (resultprime - result)/dsigma;
    end
    
%     diff = diff + eye(M)*realmin;
end
