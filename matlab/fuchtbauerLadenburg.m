function scale = fuchtbauerLadenburg(tau, sigmaFun)
%FUCHTBAUERLADENBURG Use Füchtbauer-Ladenburg equation to scale emission sigma
%   Returns the factor the current sigmaemi_* needs to be scaled by

% Constants
c = 3e8;% m/s

% Interval of integration (must cover the entire function)
minlambda = 240e-9; % m
maxlambda = 740e-9; % m

ncore = @refractionIndexPMMA;

integrand = @(ll) ncore(ll).^2.*sigmaFun(ll)./ll.^4;

scale = 1/(8*pi*c*tau*quadgk(integrand, minlambda, maxlambda));

end