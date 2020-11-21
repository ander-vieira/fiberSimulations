function sigmaabs = sigmaabs_AC56(lambda)
%SIGMAABS_AC46 Summary of this function goes here
%   Detailed explanation goes here

referenceLambda = 379e-9;
referenceSigma = 6.9e-23;

rawdata = [               % primera columna lambda (nm), segunda columna sigma emisión sin normalizar (a.u.)    
    230        0
    235        0
    237.18  0
    249.38  0
  221.67            0
       380.23            0
          384   1.3381e-22
       387.04   3.5184e-22
       390.07   7.5779e-22
       396.88   3.0041e-21
       401.04   4.8175e-21
       406.72   6.6849e-21
       411.26   8.2547e-21
          418   1.0068e-20
          427   1.2044e-20
       435.86   1.3722e-20
       439.26   1.4263e-20
       441.53   1.4452e-20
       444.94   1.4561e-20
        449.1   1.4588e-20
       452.13   1.4506e-20
       456.67   1.4452e-20
        459.7   1.4506e-20
       462.72   1.4669e-20
          465   1.4777e-20
       467.27   1.4777e-20
       469.16   1.4723e-20
       471.81   1.4425e-20
       475.21   1.3316e-20
       479.75   1.0744e-20
       484.67   7.3073e-21
       488.46   4.4385e-21
       491.86   2.1651e-21
       495.65   1.0014e-21
        498.3   5.9542e-22
       499.81   4.0597e-22
       502.08   2.1651e-22
          506   1.0035e-22
       516.46            0
       699.62            0
    752               0
    753               0
    754               0
    755               0
    ];


rawdata(:,2) = rawdata(:,2)*0.96113;
% rawdata = csvread('sigmaabs_AC56.csv');
rawlambdas = rawdata(:, 1)*1e-9;
rawvalues = rawdata(:, 2);

% rawvalues = rawvalues*referenceSigma/spline(rawlambdas, rawvalues, referenceLambda);

sigmaabs = spline(rawlambdas, rawvalues, lambda);

sigmaabs = sigmaabs.*(sigmaabs>=0).*(lambda>=min(rawlambdas)).*(lambda<=max(rawlambdas));

end

