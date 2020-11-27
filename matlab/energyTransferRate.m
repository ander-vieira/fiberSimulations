function rate = energyTransferRate(donorDopant, acceptorDopant, acceptorN)
%ENERGYTRANSFERRATE Summary of this function goes here
%   Detailed explanation goes here

ll = (240:1:740)*1e-9;

[tau, ~, donorEmiFun] = getDyeDopantAttributes(donorDopant);
[~, acceptorAbsFun, ~] = getDyeDopantAttributes(acceptorDopant);

kappa2 = 2/3;
quantumYield = 1;

donorEmi = donorEmiFun(ll);
donorWNSP = donorEmi/sum(donorEmi);
acceptorAbs = acceptorAbsFun(ll);

nPMMA = refractionIndexPMMA(ll);

integral = sum(donorWNSP.*acceptorAbs.*ll.^4./nPMMA.^4);

rate = acceptorN^2*9*log(10)*kappa2*quantumYield*integral/(1280*pi^5)/tau;

end

