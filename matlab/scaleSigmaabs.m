function scale = scaleSigmaabs(sigmaabs, tau)

c = 3e8; % m/s
dielectricConstant = 8.85e-12; % C^2*s^2*kg^-1*m^-3

electronCharge = 1.6e-19; % C
electronMass = 9.1e-31; % kg

% One method: using average wavelength values
sigmaabsW = @(w) sigmaabs(2*pi*c./w);

avgW = quadgk(@(w) w.*sigmaabsW(w), 0, 1e16)/quadgk(sigmaabsW, 0, 1e16);
avgLambda = 2*pi*c/avgW;

n = refractionIndexPMMA(avgLambda);

gammaBase = n*electronCharge^2*avgW^2/(6*pi*electronMass*c^3*dielectricConstant);
dopantStrength = 1/(gammaBase*tau);

scale = pi*electronCharge^2*dopantStrength/(2*electronMass*c*dielectricConstant*n*quadgk(sigmaabsW, 0, 1e16));

% Another method: put all the lambda-dependent values into the integral
% (see Fuchtbauer-Ladenburg)
integrandLambda = @(lambda) sigmaabs(lambda).*refractionIndexPMMA(lambda).^2./lambda.^4; % m^2

scale = 3/(8*pi*c*tau*quadgk(integrandLambda, 0, 740e-9));

end