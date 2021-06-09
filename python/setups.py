import math;

import vector3;
import params;
import constraints;
import interphases;
import components;

def noClad(dopant, N, diameter, lightL, darkL, ll):
    generator = components.SolarGenerator(params.solarIrradiance, ll, vector3.new(0, -1, 0), vector3.new(2*diameter, 0, 0), vector3.new(0, 0, lightL), origin=vector3.new(-2*diameter/2, diameter+1e-5, 0));

    cint1 = constraints.Interval(0, lightL+darkL);
    cint2 = constraints.Cylinder(diameter/2);
    cint3 = constraints.join(cint1, cint2);
    #cint4 = constraints.Interval(-math.inf, 0, axis=1);
    collector = components.Collector(interphases.PlaneInterphase(vector3.new(0, 0, 1), origin=vector3.new(0, 0, lightL+darkL)), len(ll), cint=cint2);
    cmp1 = components.Refractor(interphases.CylinderInterphase(diameter/2), "air", "PMMA", ll, cint=cint1);
    cmp2 = components.Attenuation("PMMA", ll, cint=cint3);
    cmp3 = components.DyeDopant("Rh6G", 1.5e22, ll, cint=cint3);
    #cmp4 = components.Mirror(interphases.PlaneInterphase(vector3.new(0, 1, 0), origin=vector3.new(0, -diameter/2-1e-6, 0)));
    #cmp4 = components.Mirror(interphases.ParabolaInterphase(3*diameter/4, origin=vector3.new(0, -diameter/2-1e-6, 0)));
    #cmp4 = components.Mirror(interphases.CylinderInterphase(diameter), cint=cint4);
    #cmp5 = components.Mirror(interphases.PlaneInterphase(vector3.new(0, 0, 1)), cint=cint2);
    #cmpList = [collector, cmp1, cmp2, cmp3, cmp4, cmp5];
    cmpList = [collector, cmp1, cmp2, cmp3];

    return generator, collector, cmpList;

def coreClad(dopant, N, diameter, q, lightL, darkL, ll):
    generator = components.SolarGenerator(params.solarIrradiance, ll, vector3.new(0, -1, 0), vector3.new(2*diameter, 0, 0), vector3.new(0, 0, lightL), origin=vector3.new(-2*diameter/2, diameter+1e-5, 0));

    cint1 = constraints.Interval(0, lightL+darkL);
    cint2 = constraints.Cylinder(diameter/2);
    cint3 = constraints.Cylinder(diameter*q/2);
    cint4 = lambda pos: cint1(pos) and cint3(pos);
    cint5 = lambda pos: cint1(pos) and cint2(pos) and not cint3(pos);
    collector = components.Collector(interphases.PlaneInterphase(vector3.new(0, 0, 1), origin=vector3.new(0, 0, lightL+darkL)), len(ll), cint=cint2);
    cmp1 = components.Refractor(interphases.CylinderInterphase(diameter/2), "air", "clad", ll, cint=cint1);
    cmp2 = components.Refractor(interphases.CylinderInterphase(diameter*q/2), "clad", "PMMA", ll, cint=cint1);
    cmp3 = components.Attenuation("PMMA", ll, cint=cint4);
    cmp4 = components.DyeDopant(dopant, N, ll, cint=cint4);
    cmp5 = components.Attenuation("clad", ll, cint=cint5);
    cmpList = [collector, cmp1, cmp2, cmp3, cmp4, cmp5];

    return generator, collector, cmpList;

def dopantInClad(dopant, N, diameter, q, lightL, darkL, ll):
    generator = components.SolarGenerator(params.solarIrradiance, ll, vector3.new(0, -1, 0), vector3.new(2*diameter, 0, 0), vector3.new(0, 0, lightL), origin=vector3.new(-2*diameter/2, diameter+1e-5, 0));

    cint1 = constraints.Interval(0, lightL+darkL);
    cint2 = constraints.Cylinder(diameter/2);
    cint3 = constraints.Cylinder(diameter*q/2);
    cint4 = lambda pos: cint1(pos) and cint3(pos);
    cint5 = lambda pos: cint1(pos) and cint2(pos) and not cint3(pos);
    collector = components.Collector(interphases.PlaneInterphase(vector3.new(0, 0, 1), origin=vector3.new(0, 0, lightL+darkL)), len(ll), cint=cint2);
    cmp1 = components.Refractor(interphases.CylinderInterphase(diameter/2), "air", "PMMA", ll, cint=cint1);
    cmp2 = components.Refractor(interphases.CylinderInterphase(diameter*q/2), "PMMA", "PMMA", ll, cint=cint1);
    cmp3 = components.Attenuation("PMMA", ll, cint=cint4);
    cmp4 = components.DyeDopant(dopant, N, ll, cint=cint5);
    cmp5 = components.Attenuation("PMMA", ll, cint=cint5);
    cmpList = [collector, cmp1, cmp2, cmp3, cmp4, cmp5];

    return generator, collector, cmpList;
