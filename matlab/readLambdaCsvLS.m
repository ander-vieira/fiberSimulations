function values = readLambdaCsvLS(csvFile, lambdas, referenceLambda, referenceValue)
%READLAMBDACSVLS Summary of this function goes here
%   Detailed explanation goes here

% Read data from the CSV file
rawData = csvread(csvFile);
rawLambdas = rawData(:, 1)*1e-9;
rawValues = rawData(:, 2);

% Scale data using the reference values
rawValues = rawValues*referenceValue/interpolateLS(rawLambdas, rawValues, referenceLambda);

values = zeros(1, length(lambdas));
for i = 1:length(lambdas)
    values(i) = interpolateLS(rawLambdas, rawValues, lambdas(i));
end

end

function value = interpolateLS(rawLambdas, rawValues, lambda)
    windowSize = 3;

    if lambda <= min(rawLambdas) || lambda > max(rawLambdas)
        value = 0;
    else
        lesserValues = find(rawLambdas<lambda);
        lesserIndex  = lesserValues(end);
        
        if lesserIndex - windowSize > 0 && lesserIndex + windowSize < length(rawLambdas)
            u = (lambda-rawLambdas(lesserIndex))/(rawLambdas(lesserIndex+1)-rawLambdas(lesserIndex));
            
            lambdaLesser = rawLambdas(lesserIndex-windowSize)*(1-u)+rawLambdas(lesserIndex-windowSize+1)*u;
            valueLesser = rawValues(lesserIndex-windowSize)*(1-u)+rawValues(lesserIndex-windowSize+1)*u;
            lambdaGreater = rawLambdas(lesserIndex+windowSize)*(1-u)+rawLambdas(lesserIndex+windowSize+1)*u;
            valueGreater = rawValues(lesserIndex+windowSize)*(1-u)+rawValues(lesserIndex+windowSize+1)*u;
            
            lambdaSample = [lambdaLesser ; rawLambdas((lesserIndex-windowSize+2):(lesserIndex+windowSize-1)) ; lambdaGreater];
            valueSample = [valueLesser ; rawValues((lesserIndex-windowSize+2):(lesserIndex+windowSize-1)) ; valueGreater];
            
            A = [ones(size(lambdaSample)) lambdaSample];
            b = valueSample;
            
            v = A\b;
            
            value = v(2)*lambda+v(1);
        else
            value = (rawValues(lesserIndex)+rawValues(lesserIndex+1))/2;
        end
    end
end