function [curve] = lambdaLogNormal(ll, lambdaIn, sigmaabs, sigmaemi)
%LAMBDALOGNORMAL Summary of this function goes here
%   Detailed explanation goes here

% Calculate expected value and deviation of average losses (in m^-1,
% scaled down by a factor of h*c)
lossAverage = sum((1./ll-1./ll').*sigmaabs.*sigmaemi', 'all')/sum(sigmaabs)/sum(sigmaemi); % m^-1
lossMoment2 = sum((1./ll-1./ll').^2.*sigmaabs.*sigmaemi', 'all')/sum(sigmaabs)/sum(sigmaemi); % m^-2

% Get parameters for log-normal distribution
mu = log(lossAverage^2/sqrt(lossMoment2));
sigma = sqrt(log(lossMoment2/lossAverage^2));

% Generate curve using log-normal distribution
curve = exp(-((log(abs(1./lambdaIn-1./ll))-mu)/sigma).^2/2)./(ll.^2./lambdaIn-ll)/sqrt(2*pi)/sigma.*(ll>lambdaIn);

end

