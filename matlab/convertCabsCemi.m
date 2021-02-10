function [cValsemi] = convertCabsCemi(lambdaVals, cValsabs, deltaWVals, tau, T)
%CONVERTCABSCEMI Convert gaussian weights from absorption to emission
%   An absorption spectrum can be converted into a sum of gaussian curves
%   using gaussianSigmasMC. The emission spectrum 

h = 6.63e-34; %  j·s
c = 3e8; % m/s
kB = 1.38e-23; % J/K

cValsemi = cValsabs.*exp(h*c/kB/T./lambdaVals);

sigmaemi = generateGaussianSigma(lambdaVals, cValsemi, deltaWVals);

cValsemi = cValsemi * fuchtbauerLadenburg(tau, sigmaemi);

end

