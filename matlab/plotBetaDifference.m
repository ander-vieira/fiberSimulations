close all;
clear;

ll = (240:5:740)*1e-9;
numll = length(ll);

ncore = refractionIndexPMMA(ll);

betaOld = (1-1./ncore)/2;

betaIntegral = zeros(1, numll);
betaEasy = zeros(1, numll);
betaRaytracing = zeros(1, numll);

for i=1:numll
    betaIntegral(i) = calculateBetaIntegral(ll(i));
    betaEasy(i) = calculateBetaEasy(ll(i));
    betaRaytracing(i) = calculateBetaRaytracing(1e-3, 0.001, ll(i));
end

figure(1);
plot(ll*1e9, betaOld);
hold on;
plot(ll*1e9, betaIntegral);
plot(ll*1e9, betaEasy);
plot(ll*1e9,betaRaytracing);
title('Beta (usable fraction of spontaneous emission)');
legend('Old method (volume of cone)', 'New method, via integral', 'New method, via raytracing 1', 'New method, via raytracing 2');
xlabel('\lambda (nm)');