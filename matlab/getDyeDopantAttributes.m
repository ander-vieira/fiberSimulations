function [tauRad, sigmaabs, sigmaemi, tauNR, quenchingA, quenchingB] = getDyeDopantAttributes(dopant)

if strcmp(dopant, 'RhB')
    tauRad = 4.8e-9;
    tauNR = 1; % Very large value -> quantum yield is basically 1
    sigmaabs = readLambdaCsv("../csv/sigmaabs_RhB.csv", 559.29e-9, 3.37e-20);
    sigmaemi = readLambdaCsv("../csv/sigmaemi_RhB.csv", 571.08e-9, 2.4973e-20);
    quenchingA = 0; % Not known -> no quenching happens
    quenchingB = -3;
    
elseif strcmp(dopant, 'C1')
    tauRad = 3.4e-9;
    tauNR = 1; % Very large value -> quantum yield is basically 1
    sigmaabs = readLambdaCsv("../csv/sigmaabs_C1.csv", 341.01e-9, 8.3465e-21);
    sigmaemi = readLambdaCsv("../csv/sigmaemi_C1.csv", 404.84e-9, 8.3713e-21);
    quenchingA = 0; % Not known -> no quenching happens
    quenchingB = -3;
    
elseif strcmp(dopant, 'C6')
    tauRad = 6.8e-9; % Need to check again
    tauNR = 2.41e-8; % Calculated from quantum yield of 0.78
%     sigmaabs = readLambdaCsv("../csv/sigmaabs_C6.csv", 459e-9, 7.290817e-21/2);
%     sigmaemi = readLambdaCsv("../csv/sigmaemi_C6.csv", 501e-9, 1.183142e-20*1.01768/2);
    sigmaabs = generateGaussianSigma([425.58e-9;449.42e-9;473.00e-9], [5.34569e-21;1.008615e-20;4.559433e-21], [4.36275e14;2.43868e14;1.19762e14]);
    sigmaemi = generateGaussianSigma([494.76e-9;514.39e-9;552.96e-9]-52e-9, [3.39061e-21;4.71290e-21;1.37855e-21], [1.3014e14;2.40312e14;3.02088e14]);
    quenchingA = 0; % Not known -> no quenching happens
    quenchingB = -3;
    
elseif strcmp(dopant, 'LumogenRed')
    tauRad = 6e-9;
    tauNR = 1; % Very large value -> quantum yield is basically 1
    sigmaabs = readLambdaCsv("../csv/sigmaabs_LumR.csv", 573.37e-9, 1.5382e-20);
    sigmaemi = readLambdaCsv("../csv/sigmaemi_LumR.csv", 600.85e-9, 1.532e-20);
    quenchingA = 7.3e9;
    quenchingB = -3.92;
    
elseif strcmp(dopant, 'LumogenOrange')
    tauRad = 6e-9;
    tauNR = 1.14e-7; % Calculated from a quantum yield of 0.95
%     sigmaabs = readLambdaCsv("../csv/sigmaabs_LumO.csv", 527e-9, 2.5290e-20);
%     sigmaemi = readLambdaCsv("../csv/sigmaemi_LumO.csv", 569e-9, 2.1047e-20);
    sigmaabs = generateGaussianSigma([411.93e-9;452.19e-9;487.83e-9;518.51e-9;525.17e-9], [3.80555e-21;4.18278e-21;1.08276e-20;4.03081e-21;1.93256e-20], [3.98957e14;1.04652e14;1.53683e14;4.57131e14;8.45222e13]);
    sigmaemi = generateGaussianSigma([566.04e-9;572.77e-9;620.67e-9]-20e-9, [1.23702e-20;9.58196e-21;4.79994e-21], [1.50599e14;7.84788e13;1.3689e14]);
    quenchingA = 0;
    quenchingB = -3;
    
else % Default is Rh6G
    tauRad = 4.8e-9;
    tauNR = 1; % Very large value -> quantum yield is basically 1
    sigmaabs = readLambdaCsv("../csv/sigmaabs_Rh6G.csv", 530.58e-9, 4.3298e-20);
    sigmaemi = readLambdaCsv("../csv/sigmaemi_Rh6G.csv", 544.02e-9, 2.0375e-20);
%     sigmaabs = generateGaussianSigma([508.84e-9;530.03e-9], [1.62791e-20;2.74366e-20], [4.01125e14;1.41498e14]);
%     sigmaemi = generateGaussianSigma([541.52e-9;577.56e-9], [1.67023e-20;5.16708e-21], [1.65728e14;2.47508e14]);
    quenchingA = 7e9;
    quenchingB = -2.4;
end

end