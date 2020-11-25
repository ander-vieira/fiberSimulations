function sigmaemi = sigmaemi_AC46(lambda)
%SIGMAABS_AC46 Summary of this function goes here
%   Detailed explanation goes here

referenceLambda = 614e-9;
referenceSigma = 4.64e-24;

rawdata = csvread('../csv/sigmaemi_AC46.csv');
rawlambdas = rawdata(:, 1)*1e-9;
rawvalues = rawdata(:, 2);

rawvalues = rawvalues*referenceSigma/spline(rawlambdas, rawvalues, referenceLambda);

sigmaemi = spline(rawlambdas, rawvalues, lambda);

sigmaemi = sigmaemi.*(sigmaemi>=0).*(lambda>=min(rawlambdas)).*(lambda<=max(rawlambdas));

end

