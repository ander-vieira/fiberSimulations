function scale = fuchtbauerLadenburg(dopant)

c = 3e8;% m/s

[tau, ~, sigmaemiFun] = getDyeDopantAttributes(dopant);

minlambda = 240e-9; % m
maxlambda = 740e-9; % m

ncore = @refractionIndexPMMA;

integrand = @(ll) ncore(ll).^2.*sigmaemiFun(ll)./ll.^4;

scale = 1/(8*pi*c*tau*quadgk(integrand, minlambda, maxlambda));

end