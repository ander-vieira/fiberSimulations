function result = readLambdaCsvLinear(csvFile, referenceLambda, referenceValue)
%READLAMBDACSVLINEAR Summary of this function goes here
%   Detailed explanation goes here

% Read data from the CSV file
rawData = csvread(csvFile);
rawLambdas = rawData(:, 1)*1e-9;
rawValues = rawData(:, 2);

function values = resultFun(lambdas)
    if(length(lambdas) > 1)
        values = zeros(1, length(lambdas));
        
        for i=1:length(lambdas)
            values(i) = resultFun(lambdas(i));
        end
    else
        indexLeft = nnz(rawLambdas < lambdas);
        
        if indexLeft == 0 || indexLeft == length(rawLambdas)
            values = 0;
        else
            lambdaLeft = rawLambdas(indexLeft);
            lambdaRight = rawLambdas(indexLeft+1);
            valueLeft = rawValues(indexLeft);
            valueRight = rawValues(indexLeft+1);

            u = (lambdas-lambdaLeft)/(lambdaRight-lambdaLeft);
            values = valueLeft*(1-u)+valueRight*u;
        end
    end
end

if resultFun(referenceLambda) ~= 0
    rawValues = rawValues*referenceValue/resultFun(referenceLambda);
end

result = @resultFun;

end

