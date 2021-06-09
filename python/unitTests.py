import vector3;
import basics;
import constraints;
import interphases;
import components;
import sims;

diameter = 1e-3;
L = 0.1;

ll = sims.generateLambdas(440e-9, 740e-9, 101);

cmp1 = components.Refractor(interphases.CylinderInterphase(diameter/2), "air", "PMMA", ll);
cmp1.setConstraint(constraints.Interval(0, L));
cmp1.process(basics.Ray(vector3.new(0, diameter/2, -L), vector3.new(0, 0.8, 0.6), 4, 1), 1e-5);

collector = components.Collector(interphases.PlaneInterphase(vector3.new(0, 0, 1)), len(ll));
collector.setConstraint(constraints.Cylinder(diameter/2));
print(collector.process(basics.Ray(vector3.new(0, 0, L), vector3.new(0, -0.8, 0.6), 4, 0.5), 1e-5));
print(sum(collector.finalPower));
