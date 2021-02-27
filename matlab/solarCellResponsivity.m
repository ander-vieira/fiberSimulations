function solarCellR = solarCellResponsivity(ll)
%SOLARCELLRESPONSIVITY Summary of this function goes here
%   Detailed explanation goes here

solarCellR = readLambdaCsv("../csv/solarCellResponsivity.csv", ll, 930.45e-9, 0.6478);

end

