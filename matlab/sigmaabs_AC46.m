function sigmaabs = sigmaabs_AC46(lambda)
%SIGMAABS_AC46 Summary of this function goes here
%   Detailed explanation goes here

referenceLambda = 379e-9;
referenceSigma = 3.45e-23;

rawdata = csvread('../csv/sigmaabs_AC46.csv');
rawlambdas = rawdata(:, 1)*1e-9;
rawvalues = rawdata(:, 2);

rawvalues = rawvalues*referenceSigma/spline(rawlambdas, rawvalues, referenceLambda);

sigmaabs = spline(rawlambdas, rawvalues, lambda);

sigmaabs = sigmaabs.*(sigmaabs>=0).*(lambda>=min(rawlambdas)).*(lambda<=max(rawlambdas));

end

