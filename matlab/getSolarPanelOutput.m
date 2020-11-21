clear;

minlambda = 240e-9;
dlambda = 1e-9;
maxlambda = 1100e-9;
ll = minlambda:dlambda:maxlambda;

weight = 20; % kg
solarType = "space";

isol = solarIrradianceSpline(ll, solarType);
panelDensity = 4.12; % kg/m^2
panelArea = weight / panelDensity; % m^2
Ppanel_in = isol*dlambda*panelArea; % W
Ppanel_out = solarCellConversion2(ll, Ppanel_in, panelArea); % W
solarConstant = 1361; % W/m^2
Psilicon_out = solarConstant*panelArea*0.169; % W
Ptriple_out = solarConstant*panelArea*0.33; % W

diameter = 0.1e-3; % m
concDensity = 1.5*diameter*1.18e3; % kg/m^2, for thickness=diameter, +20% for other components
concArea = weight / concDensity; % m^2
L = 0.2; % m
W = concArea / L; % m
fiberN = W/diameter;
[~, Pfiber_out] = combinedDopantIterative(["C6"], [7e23], ["AC56"], [1e22], solarType, diameter, L); 
Pconc_out = Pfiber_out*fiberN*2*1.6;

fprintf("Concentrator density: %g kg/m^2\n", concDensity);
fprintf("Simulated solar cell output power: %g W (area %g m^2)\n", Ppanel_out, panelArea);
fprintf("Silicon solar cell output power: %g W (area %g m^2)\n", Psilicon_out, panelArea);
fprintf("Triple junction solar cell output power: %g W (area %g m^2)\n", Ptriple_out, panelArea);
fprintf("Concentrator output power: %g W (area %g m^2)\n", Pconc_out, concArea);
fprintf("Fraction of maximum: %.5f\n", Pconc_out/Ptriple_out);