function [beta, Kz] = geometricalParamsB(ncore)
%GEOMETRICALPARAMSB Summary of this function goes here
%   Detailed explanation goes here

beta = (ncore-1)/ncore/2;
Kz = 2*ncore/(ncore+1);

end

