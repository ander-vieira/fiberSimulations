function alpha = attenuationPMMA(ll)
%ATTENUATIONPMMA Summary of this function goes here
%   Detailed explanation goes here

alphaFun = readLambdaCsv("../csv/attenuationPMMA.csv", 500e-9, 0.01842); % Scaled for 0.08 dB/m at green (500 nm)

alpha = alphaFun(ll);

end

