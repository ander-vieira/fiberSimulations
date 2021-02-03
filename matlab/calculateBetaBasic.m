function beta = calculateBetaBasic(ll)
%CALCULATEBETABASIC Summary of this function goes here
%   Detailed explanation goes here

ncore = refractionIndexPMMA(ll);

beta = (ncore - 1)./(2*ncore);

end