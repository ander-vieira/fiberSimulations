function irr = solarIrradianceSpline(lambda, solarType, rawdata)

if nargin < 2
    solarType = "surface";
end
if(nargin < 3 || isempty(rawdata))
    % Read data from the csv file
	rawdata = csvread('../csv/solarIrradiance.csv');
end
rawlambdas = rawdata(:, 1)*1e-9;

if strcmp(solarType, "space")
    rawvalues = rawdata(:, 2);
    rawvalues = rawvalues * 1361 / sum(rawvalues);
elseif strcmp(solarType, "spacemars")
    rawvalues = rawdata(:, 2);
    rawvalues = rawvalues * 587 / sum(rawvalues);
elseif strcmp(solarType, "spacevenus")
    rawvalues = rawdata(:, 2);
    rawvalues = rawvalues * 2604 / sum(rawvalues);
else
    rawvalues = rawdata(:, 3);
    % Normalize total irradiance to 900 W/m^2
    rawvalues = rawvalues * 1000 / sum(rawvalues);
end

% Check for the required parameter
if nargin == 0
    figure(1);
    plot(rawlambdas*1e9, rawvalues);
    title('Spectrum of solar irradiance (source: ASTM G-173)');
    xlabel('Wavelength (nm)');
    ylabel('I_sol (W/m^2/nm)');
    disp('Returns value of standard solar irradiance in W/(m^2*nm)');
    disp('Lambda (scalar or vector parameter) is required.');
    
    return
end

%Lambda range check
minlambda = 280e-9;
maxlambda = 4000e-9;

% Low pass filter parameters
% If either is 0, there's no LPF applied
dlambda = 0.5e-9;
windowSize = 2;

% Spline interpolation N times, then averages all of them (low pass filter)
lambdaM = (-windowSize:windowSize)'*dlambda.*ones(1, length(lambda))+lambda;
irrM = spline(rawlambdas, rawvalues, lambdaM);
irr = sum(irrM, 1)/(2*windowSize+1);
% irr = spline(rawlambdas, rawvalues, lambda);

% Remove negatives + stuff outside limits
irr = irr.*(irr>=0).*(lambda>=minlambda).*(lambda<=maxlambda);

% Scale result to W/m^3 from W/(m^2*nm)
irr = irr*1e9;

return
