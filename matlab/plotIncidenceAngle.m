clear;
close all;
tic;

angles = [0 10 20 30 40 50]; % º

dopant = "LumogenRed";
N = 1.508e22;
diameter = 1e-3;
lightL = 0.05;
darkL = 0;

M = length(angles);
outputPower = zeros(1, M);

for i = 1:M
    fprintf('%d/%d\n', i, M);
    
    outputPower(i) = oneDopantRaytracing(dopant, N, diameter, lightL, darkL, angles(i));
end

fprintf('Simulation time: %.1f s\n', toc());
figure(1);
plot(angles, outputPower*1e6);
xlabel('Incidence angle (º)');
ylabel('Output power (\muW)');
title('Output power over incidence angle of sunlight (with raytracing)');