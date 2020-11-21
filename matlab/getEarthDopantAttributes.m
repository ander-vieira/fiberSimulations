function [tauT, tauD, wTD, wDT, sigmaabs, sigmaemiD] = getEarthDopantAttributes(dopant)

if strcmp(dopant, 'AC56')
    tauT = 1e-5;
    tauD = 1.02e-3;
    wTD = 8e8;
    wDT = 2e8;
    sigmaabs = @sigmaabs_AC56;
    sigmaemiD = @sigmaemi_AC56;
else % Default is AC46
    tauT = 1e-5;
    tauD = 1.02e-3;
    wTD = 8e8;
    wDT = 2e8;
    sigmaabs = @sigmaabs_AC46;
    sigmaemiD = @sigmaemi_AC46;
end

end