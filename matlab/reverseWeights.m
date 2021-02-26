function [w] = reverseWeights(N, lightL, darkL, ll, sigmaabs, Pout)
%REVERSEWEIGHTS Obtain values of the WNSP from measured output power

tic;

ncore = refractionIndexPMMA(ll);
alfaPMMA = valuesalfaPMMA(ll);

beta = (ncore - 1)./(2*ncore);

% Interval used for numerical differentiation
dw = 1e-7;

% Normalize since only the shape is needed
normalizedP = Pout/sum(Pout);

% Iterate using Newton-Raphson adapted for M-sized vectors
w = normalizedP;
error = 1;
while error > 1e-5
    % Get expected power spectrum shape
    [normalizedFun, J] = solveSystem(ll, N, lightL, darkL, sigmaabs, alfaPMMA, beta, w, dw);
    
    A = [J'*J ones(length(ll),1) ; ones(1, length(ll)) 0];
    
    b = [J'*(normalizedP'-normalizedFun) ; 0];
    
    c = (A\b)';
    deltaW = c(1:end-1);
    
    % Calculate relative error in each step
    error = norm(deltaW);
    
    % Calculate new values for the WNSP
    w = w + 0.1*deltaW;
end

% Remove negative values
w(w < 0) = 0;
w = w / sum(w);
    
fprintf("Time elapsed in reverseWeights: %.2f\n", toc());
    
end

% Simulates the behaviour of the POF via the stationary method
% Returns the shape of the output power spectrum
% Also returns the numerical differential of the shape of P
% with respect to each value in WNSP for N-R purposes
function [result, diff] = solveSystem(ll, N, lightL, darkL, sigmaabs, alfaPMMA, beta, wnsp, dw)
    M = length(ll);
    diff = zeros(M);
    
    % Define C and G, which constitute a system of linear differential
    % equations
    C = beta.*wnsp./ll;

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
    result = result/sum(result);
    
    % Get the differential for each component k of WNSP
    for k = 1:M
        % Increase component k of WNSP slightly, normalize
        direction = [zeros(1,k-1) 1 zeros(1,M-k)];
        newW = (wnsp+dw*direction)/(1+dw);
        
        % Define a new C and G for the modified WNSP
        Cprime = beta.*newW./ll;
    
        for i = 1:M
            for j = 1:M
                if(i == j)
                    Gprime(i, i) = -sigmaabs(i).*N.*(1-beta(i).*newW(i))-alfaPMMA(i)+sqrt(realmin);
                else
                    Gprime(j, i) = sigmaabs(i)*N*beta(j)*newW(j)*ll(i)/ll(j);
                end
            end
        end
        
        % Solve the system of linear differential equations again
        resultprime = expm(Gprime*darkL)*(expm(Gprime*lightL)-eye(M))*(Gprime\Cprime');
        resultprime = resultprime/sum(resultprime);
        
        % Get the derivative from substracting both results
        diff(:, k) = (resultprime - result)/dw;
    end
    
    diff = diff + eye(M)*sqrt(realmin);
end