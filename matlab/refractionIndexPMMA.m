function refractionIndex = refractionIndexPMMA(ll)
%REFRACTIONINDEXPMMA Summary of this function goes here
%   Detailed explanation goes here

refractionIndex = sqrt(1./(-8226e-18./ll.^2+0.8393)+1);

end

