function plotDyeDopantSigmas(dopant, figura)

if nargin < 2
    figura = 1;
end

[~, sigmaabsFun, sigmaemiFun] = getDyeDopantAttributes(dopant);

ll = (240:1:840)*1e-9;

sigmaabs = sigmaabsFun(ll);
sigmaemi = sigmaemiFun(ll);

close(figura);
figure(figura);
plot(ll*1e9, sigmaabs);
hold on;
plot(ll*1e9, sigmaemi);
xlabel('\lambda (nm)');
ylabel('\sigma (m^2)');
title(sprintf('Sigmas for %s', dopant));

end