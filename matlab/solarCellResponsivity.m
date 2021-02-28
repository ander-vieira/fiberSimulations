function solarCellR = solarCellResponsivity(ll)
%SOLARCELLRESPONSIVITY Summary of this function goes here
%   Detailed explanation goes here

solarCellFun = readLambdaCsv("../csv/solarCellResponsivity.csv", 930.45e-9, 0.6478);

solarCellR = solarCellFun(ll);

end

