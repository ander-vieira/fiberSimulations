function [tau, sigmaabs, sigmaemi, quenchingA, quenchingB] = getDyeDopantAttributes(dopant)

if strcmp(dopant, 'RhB')
    tau = 4.8e-9;
    sigmaabs = @sigmaabs_RhB;
    sigmaemi = @sigmaemi_RhB;
    quenchingA = 0; % Not known -> no quenching happens
    quenchingB = -3;
elseif strcmp(dopant, 'C1')
    tau = 3.4e-9;
    sigmaabs = @(lambdas) readLambdaCsv("../csv/sigmaabs_C1.csv", lambdas, 341.01e-9, 8.3465e-21);
    sigmaemi = @(lambdas) readLambdaCsv("../csv/sigmaemi_C1.csv", lambdas, 404.84e-9, 8.3713e-21);
    quenchingA = 0; % Not known -> no quenching happens
    quenchingB = -3;
elseif strcmp(dopant, 'C6')
    tau = 3.4e-9; % Need to check again
    sigmaabs = @sigmaabs_C6;
    sigmaemi = @sigmaemi_C6;
    quenchingA = 0; % Not known -> no quenching happens
    quenchingB = -3;
elseif strcmp(dopant, 'LumogenRed')
    tau = 6e-9;
    sigmaabs = @sigmaabs_LumogenRed;
    sigmaemi = @sigmaemi_LumogenRed;
    quenchingA = 7.3e9;
    quenchingB = -3.92;
elseif strcmp(dopant, 'LumogenOrange')
    tau = 6e-9;
    sigmaabs = @(lambdas) readLambdaCsv("../csv/sigmaabs_LumO.csv", lambdas, 527e-9, 2.5351e-20);
    sigmaemi = @(lambdas) readLambdaCsv("../csv/sigmaemi_LumO.csv", lambdas, 569e-9, 2.5168e-20);
    quenchingA = 0;
    quenchingB = -3;
else % Default is Rh6G
    tau = 4.8e-9;
    sigmaabs = @sigmaabs_Rh6G;
    sigmaemi = @sigmaemi_Rh6G;
    quenchingA = 7e9;
    quenchingB = -2.4;
end

end