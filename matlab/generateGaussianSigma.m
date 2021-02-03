function sigma = generateGaussianSigma(lambdaVals, cVals, deltaWVals)
%GENERATEGAUSSIANSIGMAS Summary of this function goes here
%   Detailed explanation goes here

h = 6.63e-34; % J·s
c = 3e8; % m/S

wVals = 2*pi*c./lambdaVals;

sigmaW = @(ww) sum(cVals.*exp(-(ww-wVals).^2./deltaWVals.^2));
sigma = @(ll) sigmaW(2*pi*c./ll);

end

