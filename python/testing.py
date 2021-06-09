import numpy;

import vector3;
import basics;
import constraints;
import interphases;
import components;
import sims;

diameter = 1e-3;
q = 0.98;
L = 0.1;

minLL = 440e-9;
maxLL = 740e-9;
numLL = 101;
ll = numpy.linspace(minLL, maxLL, numLL);
dLL = ll[1]-ll[0];

cmp1 = components.Refractor(interphases.CylinderInterphase(diameter/2), "air", "PMMA", ll);
cmp1.setConstraint(constraints.Interval(0, L));
cmp1.process(basics.Ray(vector3.new(0, diameter/2, -L), vector3.new(0, 0.8, 0.6), 4, 1), 1e-5);

collector = components.Collector(interphases.PlaneInterphase(vector3.new(0, 0, 1)), numLL);
collector.setConstraint(constraints.Cylinder(diameter/2));
print(collector.process(basics.Ray(vector3.new(0, 0, L), vector3.new(0, -0.8, 0.6), 4, 0.5), 1e-5));
print(sum(collector.finalPower));

#sims.dyeRaytracing("Rh6G", 1.5e22, 3e-3, 0.98, 0.1);
sims.dyeRaytracingInClad("Rh6G", 3e22, 1e-3, 0.97, 0.1);
sims.dyeRaytracingInClad("RhB", 3e22, 1e-3, 0.97, 0.1);
sims.dyeRaytracingInClad("LumR", 3e22, 1e-3, 0.97, 0.1);
