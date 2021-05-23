function [tauT, tauD, wTD, wDT, sigmaabs, sigmaemiD] = getEarthDopantAttributes(dopant)

if strcmp(dopant, 'AC56')
    tauT = 1e-5;
    tauD = 1.02e-3;
    wTD = 8e8;
    wDT = 2e8;
    sigmaabs = readLambdaCsv("../csv/sigmaabs_AC56.csv", 379e-9, 1e-20);
    sigmaemiD = readLambdaCsv("../csv/sigmaemi_AC56.csv", 614e-9, 1e-21);
else % Default is AC46
%     tauT = 1e-5;
%     tauD = 1.02e-3;
    tauT = 1.33e-4;
    tauD = 5.33e-4;
    wTD = 8e8;
    wDT = 2e8;
    sigmaabs = readLambdaCsv("../csv/sigmaabs_AC46.csv", 379e-9, 1e-20);
    sigmaemiD = readLambdaCsv("../csv/sigmaemi_AC46.csv", 614e-9, 1e-21);
end

end