import math;
import time;
import numpy;
from matplotlib import pyplot as plt;

import setups;

def lambdaRange(minLambda, maxLambda, numLL = 61):
    return numpy.linspace(minLambda, maxLambda, numLL);

def runRay(ray, cmpList):
    while ray.loopOn:
        ds = math.inf;
        firstCmp = None;
        for cmp in cmpList:
            ds2 = cmp.intersectionPoint(ray);
            if ds2 < ds:
                ds = ds2;
                firstCmp = cmp;

        if ds < math.inf:
            ray.move(ds);
        else:
            ray.loopOn = False;

        if firstCmp is not None:
            firstCmp.process(ray, ds);
            ray.move(1e-8);

def simulate(generator, cmpList, M=200000):
    for i in range(M):
        ray = generator.generateRay(1/M);

        runRay(ray, cmpList);

        if i % (M/100) == 0:
            print("%d%%" % (i/(M/100)), end="\r");

    print("");

class resultPrinter:
    def __init__(self, collector, ll):
        self.collector = collector;
        self.ll = list(map(lambda l: l*1e9, ll));
        self.dLL = ll[1]-ll[0];
        self.startTimer();

    def startTimer(self):
        self.tic = time.time();

    def printResults(self):
        print("Elapsed time: %.2f s" % (time.time()-self.tic));
        print("Output power of fiber: %f uW" % (self.collector.getTotalPower()*1e6));

    def plotPower(self, figure=1):
        plt.figure(figure);
        plt.plot(self.ll, list(map(lambda k: 1e-3/self.dLL*self.collector.getPower(k), range(len(self.ll)))));
        plt.show();

def dyeRaytracingNoClad(dopant, N, diameter, q, lightL, darkL=0):

    ll = lambdaRange(440e-9, 740e-9, 61);

    generator, collector, cmpList = setups.noClad(dopant, N, diameter, lightL, darkL, ll);

    printer = resultPrinter(collector, ll);

    simulate(generator, cmpList, M=200000);

    printer.printResults();
    printer.plotPower(1);

def dyeRaytracing(dopant, N, diameter, q, lightL, darkL=0):

    ll = lambdaRange(440e-9, 740e-9, 61);

    generator, collector, cmpList = setups.coreClad(dopant, N, diameter, q, lightL, darkL, ll);

    printer = resultPrinter(collector, ll);

    simulate(generator, cmpList, M=200000);

    printer.printResults();
    printer.plotPower(1);

def dyeRaytracingInClad(dopant, N, diameter, q, lightL, darkL=0):

    ll = lambdaRange(440e-9, 740e-9, 61);

    generator, collector, cmpList = setups.dopantInClad(dopant, N, diameter, q, lightL, darkL, ll);

    printer = resultPrinter(collector, ll);

    simulate(generator, cmpList, M=200000);

    printer.printResults();
    printer.plotPower(1);
