function [error, bestLambda, bestC] = gaussianSigmasMC(sigmaabsFun, numVals)
% GAUSSIANSIGMASMC Get the best least-squares approximation to a sigma
% This function tries to approximate the given sigma function via a series
% of gaussian functions weighted with coefficients. If the positions of
% the gaussian curves is known, least squares approximation can be used to
% determine the coefficients.
% 
% Since knowing the wavelengths for the gaussian curves is non trivial,
% this function uses a Monte Carlo approach to get the best fitting ones:
% it generates a lot of random wavelength combinations and calculates the
% error each one yields.
% 
% Parameters:
% sigmaabsFun: a function that yields the (properly scaled) sigma to model
% numVals: how many gaussian curves should be used

if nargout == 0
    tic;
end

c = 3e8;

ll = (240:2:740)*1e-9;
ww = 2*pi*c./ll;

% Number of sets of random values to generate
N = 100000;

% Determines the width of the gaussian curves
deltaW = 2e14;

% Sample the given sigma function
sigmaabs = sigmaabsFun(ll);

% Interval of wavelengths to generate randomly between
lambdaMin = 440e-9;
lambdaMax = 640e-9;

lambdaVals = zeros(numVals, N);
cVals = zeros(numVals, N);
error = zeros(1, N);
for i = 1:N
    % Generate the random lambda values and the corresponding frequencies
    lambdaVals(:, i) = lambdaMin + (lambdaMax-lambdaMin)*rand(numVals, 1);
    wVals = 2*pi*c./lambdaVals(:, i);
    
    % Defines a matrix with the coefficients for least-squares
    % approximation. All of them are gaussian curves
    M = exp(-(ww'-wVals').^2/deltaW^2);
    
    % Use MATLAB built-in LS approximation.
    cVals(:, i) = M\sigmaabs';
    
    % Calculate the error of the approximation
    error(i) = norm(abs(sigmaabs' - M*cVals(:, i)))/norm(sigmaabs)/sqrt(N);
end

% Get the smallest error value (negative values of c are not allowed)
minError = min(error(min(cVals) >= 0));

% Get the lambda and c values for the smallest error obtained
bestLambda = lambdaVals(:, error == minError);
bestW = 2*pi*c./bestLambda;

M = exp(-(ww'-bestW').^2/deltaW^2);
bestC = cVals(:, error == minError);

if nargout == 0
    % Print figure with given and approximated sigma
    figure(1);
    plot(ll*1e9, sigmaabs);
    hold on;
    plot(ll*1e9, M*bestC);
    title(sprintf('Error: %e', min(error)));
    xlabel('\lambda (nm)');

    % Print the values that yield the best result
    fprintf('Processing time: %.1f s\n', toc());
    fprintf('Best values:\n');
    for i = 1:numVals
        fprintf('%.2f nm: c = %g\n', bestLambda(i)*1e9, bestC(i));
    end
    fprintf('Error: %g\n', minError);
end

end