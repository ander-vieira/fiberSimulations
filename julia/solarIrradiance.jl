module SolarIrradiance;

using DelimitedFiles;
using Interpolations;

solarIrradianceData = readdlm("../csv/solarIrradiance.csv", ',');
surfaceScale = 1000/sum(solarIrradianceData[:, 3]);
spaceScale = 1361/sum(solarIrradianceData[:, 2]);

rawFunction = LinearInterpolation(solarIrradianceData[:, 1]*1e-9, solarIrradianceData[:, 3]; extrapolation_bc=Throw());
rawFunctionSpace = LinearInterpolation(solarIrradianceData[:, 1]*1e-9, solarIrradianceData[:, 2]; extrapolation_bc=Throw());

surfaceValue(lambda) = rawFunction(lambda)*surfaceScale;
spaceValue(lambda) = rawFunctionSpace(lambda)*spaceScale;

end
