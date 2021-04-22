function [minLambda, maxLambda] = getLambdaRanges(dopants, dlambda)
%GETLAMBDARANGES Get the wavelength range to simulate a dopant combination
%   Each dopant's absorption/emission spectra are limited to a range of
%   wavelengths, so simulating outside that range is unnecessary. For
%   multiple dopants, the ranges of each dopant are combined. dlambda is
%   used to ensure even spacing inside the range

% Default value
if nargin < 2
    dlambda = 1e-9;
end

if length(dopants) > 1
    % Multiple dopants: call this function for each dopant, take the
    % minimum and maximum values respectively
    minLambda = 1;
    maxLambda = 0;
    
    for i=1:length(dopants)
        [minLambdaI, maxLambdaI] = getLambdaRanges(dopants(i), dlambda);
        
        minLambda = min(minLambda, minLambdaI);
        maxLambda = max(maxLambda, maxLambdaI);
    end
else
    % Values for each dopant (obtained by hand from the sigmas)
    if strcmp(dopants, "Rh6G")
        minLambda = 410e-9;
        maxLambda = 710e-9;
    elseif strcmp(dopants, "RhB")
        minLambda = 360e-9;
        maxLambda = 760e-9;
    elseif strcmp(dopants, "LumogenRed")
        minLambda = 360e-9;
        maxLambda = 750e-9;
    elseif strcmp(dopants, "C1")
        minLambda = 270e-9;
        maxLambda = 680e-9;
    elseif strcmp(dopants, "C6")
        minLambda = 320e-9;
        maxLambda = 660e-9;
    elseif strcmp(dopants, "AC46")
        minLambda = 250e-9;
        maxLambda = 750e-9;
    elseif strcmp(dopants, "AC56")
        minLambda = 250e-9;
        maxLambda = 750e-9;
    else
        minLambda = 250e-9;
        maxLambda = 750e-9;
    end
    
    % Round the values to ensure even spacing
    minLambda = dlambda*floor(minLambda/dlambda);
    maxLambda = dlambda*ceil(maxLambda/dlambda);
end

end

