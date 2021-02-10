function [tauT, tauD, wTD, wDT, sigmaabs, sigmaemiD] = getEarthDopantAttributes(dopant)

if strcmp(dopant, 'AC56')
    tauT = 1e-5;
    tauD = 1.02e-3;
    wTD = 8e8;
    wDT = 2e8;
    sigmaabs = @(lambdas) readLambdaCsv("../csv/sigmaabs_AC56.csv", lambdas, 379e-9, 6.9e-23);
    sigmaemiD = @(lambdas) readLambdaCsv("../csv/sigmaemi_AC56.csv", lambdas, 614e-9, 4.64e-24);
else % Default is AC46
    tauT = 1e-5;
    tauD = 1.02e-3;
    wTD = 8e8;
    wDT = 2e8;
    sigmaabs = @(lambdas) readLambdaCsv("../csv/sigmaabs_AC46.csv", lambdas, 379e-9, 3.45e-23);
    sigmaemiD = @(lambdas) readLambdaCsv("../csv/sigmaemi_AC46.csv", lambdas, 614e-9, 4.64e-24);
end

end