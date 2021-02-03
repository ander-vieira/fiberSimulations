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
%     sigmaabs = @sigmaabs_C6;
%     sigmaemi = @sigmaemi_C6;
%     sigmaabs = @(lambdas) readLambdaCsv("../csv/sigmaabs_C6.csv", lambdas, 459e-9, 2.8525e-20);
%     sigmaemi = @(lambdas) readLambdaCsv("../csv/sigmaemi_C6.csv", lambdas+50e-9, 501e-9, 1.8221e-20);
    sigmaabs = generateGaussianSigma([425.58e-9;449.42e-9;473.00e-9], [1.0744e-20;2.02716e-20;9.16375e-21], [4.36275e14;2.43868e14;1.19762e14]);
    sigmaemi = generateGaussianSigma([494.76e-9;514.39e-9;552.96e-9]-50e-9, [8.46871e-21;1.17714e-20;3.44319e-21], [1.3014e14;2.40312e14;3.02088e14]);
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
    sigmaabs = @(lambdas) readLambdaCsv("../csv/sigmaabs_LumO.csv", lambdas, 527e-9, 2.5290e-20);
    sigmaemi = @(lambdas) readLambdaCsv("../csv/sigmaemi_LumO.csv", lambdas, 569e-9, 2.5168e-20);
    quenchingA = 0;
    quenchingB = -3;
else % Default is Rh6G
    tau = 4.8e-9;
    sigmaabs = @sigmaabs_Rh6G;
    sigmaemi = @sigmaemi_Rh6G;
%     sigmaabs = generateGaussianSigma([508.84e-9;530.03e-9], [1.62791e-20;2.74366e-20], [4.01125e14;1.41498e14]);
%     sigmaemi = generateGaussianSigma([541.52e-9;577.56e-9], [1.67023e-20;5.16708e-21], [1.65728e14;2.47508e14]);
    quenchingA = 7e9;
    quenchingB = -2.4;
end

end