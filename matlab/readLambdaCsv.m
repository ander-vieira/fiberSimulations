function result = readLambdaCsv(csvFile, referenceLambda, referenceValue)
%READLAMBDACSV Read a function of wavelength from a CSV file
%   Reads a CSV file for lambda:value (x:y) pairs
%   and uses those values to perform spline interpolation
%   for an array of lambda values.
%   Note: lambda values in the CSV files should be in nm.
%   The function returns 0 for lambdas outside the ones in the CSV file.
%   
%   Parameters:
%   csvFile: the path to the CSV file (relative to this file or absolute)
%   lambda: array of wavelengths to get the values for (m)
%   referenceLambda: a lambda value to scale the values from the CSV
%   referenceValue: the value the function needs to have at referenceLambda

%   Example calls:
%   readLambdaCsv("../csv/sigmaabs_AC46.csv", lambdas, 379e-9, 3.45e-23);
%   readLambdaCsv("../csv/sigmaemi_AC46.csv", lambdas, 614e-9, 4.64e-24);

% Read data from the CSV file
rawData = csvread(csvFile);
rawLambdas = rawData(:, 1)*1e-9;
rawValues = rawData(:, 2);

% Scale data using the reference values
rawValues = rawValues*referenceValue/spline(rawLambdas, rawValues, referenceLambda);

% Low pass filter parameters
% If either is 0, there's no LPF applied
dlambda = 1e-9;
windowSize = 3;

% window = ones(2*windowSize+1, 1)/(2*windowSize+1);
window = exp(-(-windowSize:windowSize)'.^2/9);
window = window/sum(window);

function values = resultFun(lambdas)

% Spline interpolation N times, then averages all of them (low pass filter)
    lambdaM = (-windowSize:windowSize)'*dlambda.*ones(1, length(lambdas))+lambdas;
    valuesM = spline(rawLambdas, rawValues, lambdaM);

    % Remove negative values (!) and values outside the range of the CSV file
    valuesM = valuesM.*(valuesM>=0).*(lambdaM>=min(rawLambdas)).*(lambdaM<=max(rawLambdas));

    values = sum(valuesM.*window, 1);
end

result = @resultFun;

end

