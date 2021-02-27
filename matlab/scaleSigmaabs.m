function scale = scaleSigmaabs(tau, sigmaabsFun)
%SCALESIGMAEMI Use Füchtbauer-Ladenburg-like equation to scale absorption sigma
%   Returns the factor the current function needs to be scaled by

% Constants
c = 3e8; % m/s

% Interval of integration (must cover the entire function)
minLambda = 240e-9; % m
maxLambda = 740e-9; % m

ncore = @refractionIndexPMMA;

integrandLambda = @(lambda) sigmaabsFun(lambda).*ncore(lambda).^2./lambda.^4; % m^2

scale = 3/(8*pi*c*tau*quadgk(integrandLambda, minLambda, maxLambda));

end