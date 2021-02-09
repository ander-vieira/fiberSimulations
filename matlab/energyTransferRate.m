function rate = energyTransferRate(donorDopant, acceptorDopant, acceptorN)
%ENERGYTRANSFERRATE Calculate FRET rate between two dye dopants
%   Calculates the rate at which a dopant transfers its excited energy
%   state to another, using the two dopants' luminescence spectra and the
%   concentration of the acceptor dopant.
%   See notes for the full expression given by Forster
%
%   Parameters:
%   donorDopant: name of the donor DYE dopant
%   acceptorDopant: name of the acceptor DYE dopant
%   acceptorN: concentration of the acceptor dopant (m^-3)

ll = (240:1:740)*1e-9;

% Get the necessary parameters
[tauRad, ~, donorEmiFun] = getDyeDopantAttributes(donorDopant);
[~, acceptorAbsFun, ~] = getDyeDopantAttributes(acceptorDopant);

% Measures the degree of alignment between the dopant molecules
% A value of 2/3 is a good approximation for isotropic alignment
kappa2 = 2/3;

donorEmi = donorEmiFun(ll);
donorWNSP = donorEmi/sum(donorEmi);
acceptorAbs = acceptorAbsFun(ll);

nPMMA = refractionIndexPMMA(ll);

% Integral is done via discrete summation (could also use quadgk)
integral = sum(donorWNSP.*acceptorAbs.*ll.^4./nPMMA.^4);

% The rate is obtained from the result of the integral (s^-1)
rate = acceptorN^2*9*log(10)*kappa2*integral/(1280*pi^5)/tauRad;

end

