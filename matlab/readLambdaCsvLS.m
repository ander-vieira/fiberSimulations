function result = readLambdaCsvLS(csvFile, referenceLambda, referenceValue)
%READLAMBDACSVLS Read a function of wavelength from a CSV file
%   Ported from sigmaabs_Rh6G.m (2021-02-03), original author Felipe Jiménez
%   Reads a CSV file for lambda:value (x:y) pairs
%   and uses those values to perform windowed least squares approximation
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
rawValues = rawValues*referenceValue/interpolateLS(rawLambdas, rawValues, referenceLambda);

function values = resultFun(lambdas)
    % "Interpolate" for each individual lambda
    values = zeros(1, length(lambdas));
    for i = 1:length(lambdas)
        values(i) = interpolateLS(rawLambdas, rawValues, lambdas(i));
    end
end

result = @resultFun;

end

% This function holds the actual least squares calculation
function value = interpolateLS(rawLambdas, rawValues, lambda)
    % Should be at least 2, higher values yield smoother results
    windowSize = 3;

    if lambda <= min(rawLambdas) || lambda > max(rawLambdas)
        % Values outside the data range are 0
        value = 0;
    else
        % Get the index of the closest values in the sample data
        lesserValues = find(rawLambdas<lambda);
        lesserIndex  = lesserValues(end);
        
        if lesserIndex - windowSize > 0 && lesserIndex + windowSize < length(rawLambdas)
            % If there are enough values on both sides, do least squares
            
            % Relative distance to closest smaller value
            % Goes from 0 to 1
            % Relative distance to closest larger value is 1-u
            u = (lambda-rawLambdas(lesserIndex))/(rawLambdas(lesserIndex+1)-rawLambdas(lesserIndex));
            
            % The furthest points in the window are combined via
            % linear interpolation
            lambdaLesser = rawLambdas(lesserIndex-windowSize)*(1-u)+rawLambdas(lesserIndex-windowSize+1)*u;
            valueLesser = rawValues(lesserIndex-windowSize)*(1-u)+rawValues(lesserIndex-windowSize+1)*u;
            lambdaGreater = rawLambdas(lesserIndex+windowSize)*(1-u)+rawLambdas(lesserIndex+windowSize+1)*u;
            valueGreater = rawValues(lesserIndex+windowSize)*(1-u)+rawValues(lesserIndex+windowSize+1)*u;
            
            % The x:y (lambda:value) pairs to do least squares on
            lambdaSample = [lambdaLesser ; rawLambdas((lesserIndex-windowSize+2):(lesserIndex+windowSize-1)) ; lambdaGreater];
            valueSample = [valueLesser ; rawValues((lesserIndex-windowSize+2):(lesserIndex+windowSize-1)) ; valueGreater];
            
            % The least squares approximation is of the form p*lambda+q
            A = [ones(size(lambdaSample)) lambdaSample];
            b = valueSample;
            
            % Get the p and q coefficients (as a vector)
            v = A\b;
            
            % Get the corresponding value for the given lambda
            value = v(2)*lambda+v(1);
        else
            % If there aren't enough values, simply average the two closest
            % values
            
            value = (rawValues(lesserIndex)+rawValues(lesserIndex+1))/2;
        end
    end
end