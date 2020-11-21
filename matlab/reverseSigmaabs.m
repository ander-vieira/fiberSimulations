function [s] = reverseSigmaabs(N, L, ll, wnsp, alpha)
%REVERSEWEIGHTS Obtain values of the WNSP from measured output power

tic;

ncore = refractionIndexPMMA(ll);
alfaPMMA = valuesalfaPMMA(ll);

beta = (ncore - 1)./(2*ncore);

% Interval used for numerical differentiation
ds = 1e-25;

% Iterate using Newton-Raphson adapted for M-sized vectors
s = (alpha - alfaPMMA)/N;
error = 1;
while error > 1e-5
    [fun, J] = solveSystem(ll, N, L, alfaPMMA, beta, wnsp, s, ds);
    
    deltaS = (J\(alpha-fun)')';
    
    error = max(abs(deltaS./s));
    
    s = s + deltaS;
end

end

function [result, diff] = solveSystem(ll, N, L, alfaPMMA, beta, wnsp, sigmaabs, ds)
    M = length(ll);
    diff = zeros(M);
    
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
    
    result = -log(diag(expm(G*L))')/L;
    
    for k = 1:M
        % Increase component k of sigmaabs slightly, normalize
        direction = [zeros(1,k-1) 1 zeros(1,M-k)];
        newSigma = (wnsp+ds*direction)/(1+dw);
        
        for i = 1:M
            for j = 1:M
                if(i == j)
                    Gprime(i, i) = -newSigma(i).*N.*(1-beta(i).*wnsp(i))-alfaPMMA(i)+sqrt(realmin);
                else
                    Gprime(j, i) = newSigma(i)*N*beta(j)*wnsp(j)*ll(i)/ll(j);
                end
            end
        end
        
        resultprime = -log(diag(expm(Gprime*L))')/L;
        
        % Get the derivative from substracting both results
        diff(:, k) = (resultprime - result)/ds;
    end
end