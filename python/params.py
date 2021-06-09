import csv;
import math;
import numpy;
import scipy.interpolate;
from matplotlib import pyplot as plt;

def readCsvLL(filename, peakLL, peakValue, llCol = 0, valueCol = 1):
    reader = csv.reader(open(filename));

    llData = [];
    valueData = [];

    for row in reader:
        llData.append(float(row[llCol])*1e-9);
        valueData.append(float(row[valueCol]));

    func1 = scipy.interpolate.interp1d(llData, valueData, kind="cubic");

    minLL = min(llData);
    maxLL = max(llData);

    def func2(ll):
        dLL = 0.5e-9;
        windowSize = 2;
        funcSum = 0;
        for n in range(-windowSize, windowSize+1):
            newLL = ll+dLL*n;
            if newLL >= minLL and newLL <= maxLL and func1(newLL) >= 0:
                funcSum += func1(newLL)*1/(2*windowSize+1);

        return funcSum;

    scale = peakValue/func2(peakLL);

    return lambda ll: func2(ll)*scale;

def plotDyeDopant(dopant, minLL = 440e-9, maxLL = 740e-9):
    ll = numpy.linspace(minLL, maxLL, 201);

    sigmaabs = list(map(dyeDopants[dopant]["sigmaabs"], ll));
    sigmaemi = list(map(dyeDopants[dopant]["sigmaemi"], ll));

    plt.plot(ll, sigmaabs);
    plt.plot(ll, sigmaemi);
    plt.show();

def plotEarthDopant(dopant, minLL = 440e-9, maxLL = 740e-9):
    ll = numpy.linspace(minLL, maxLL, 201);

    sigmaabs = list(map(earthDopants[dopant]["sigmaabs"], ll));
    sigmaemi = list(map(earthDopants[dopant]["sigmaemi"], ll));

    plt.plot(ll, sigmaabs);
    plt.plot(ll, sigmaemi);
    plt.show();

def getLambdaRange(ranges):
    minLL = math.inf;
    maxLL = 0;

    for range in ranges:
        if lambdaRanges[range][0] < minLL:
            minLL = lambdaRanges[range][0];
        if lambdaRanges[range][1] > maxLL:
            maxLL = lambdaRanges[range][1];

    return minLL, maxLL;

solarIrradiance = readCsvLL("../csv/solarIrradiance.csv", 500e-9, 1.5297e9, valueCol=2);

dyeDopants = {
    "Rh6G": {
        "tauRad": 4.8e-9,
        "tauNR": 1,
        "sigmaabs": readCsvLL("../csv/sigmaabs_Rh6G.csv", 530.58e-9, 4.3298e-20),
        "sigmaemi": readCsvLL("../csv/sigmaemi_Rh6G.csv", 544.02e-9, 2.0375e-20)
    }, "RhB": {
        "tauRad": 4.8e-9,
        "tauNR": 1,
        "sigmaabs": readCsvLL("../csv/sigmaabs_RhB.csv", 559.29e-9, 3.37e-20),
        "sigmaemi": readCsvLL("../csv/sigmaemi_RhB.csv", 571.08e-9, 2.4973e-20)
    }, "LumR": {
        "tauRad": 6e-9,
        "tauNR": 1.14e-7,
        "sigmaabs": readCsvLL("../csv/sigmaabs_LumR.csv", 573.37e-9, 1.9862e-20),
        "sigmaemi": readCsvLL("../csv/sigmaemi_LumR.csv", 600.85e-9, 2.1984e-20)
    }, "LumO": {
        "tauRad": 6e-9,
        "tauNR": 1.14e-7,
        "sigmaabs": readCsvLL("../csv/sigmaabs_LumO.csv", 527.00e-9, 2.5290e-20),
        "sigmaemi": readCsvLL("../csv/sigmaemi_LumO.csv", 569.00e-9, 2.1047e-20)
    }, "LumY": {
        "tauRad": 6e-9,
        "tauNR": 7.4e-8,
        "sigmaabs": readCsvLL("../csv/sigmaabs_LumY.csv", 477.00e-9, 2.3849e-20),
        "sigmaemi": readCsvLL("../csv/sigmaemi_LumY.csv", 495.00e-9, 1.0178e-20)
    }, "C1": {
        "tauRad": 3.4e-9,
        "tauNR": 1,
        "sigmaabs": readCsvLL("../csv/sigmaabs_C1.csv", 341.01e-9, 8.3465e-21),
        "sigmaemi": readCsvLL("../csv/sigmaemi_C1.csv", 404.84e-9, 8.3713e-21)
    }, "C6": {
        "tauRad": 6.8e-9,
        "tauNR": 2.41e-8,
        "sigmaabs": readCsvLL("../csv/sigmaabs_C6.csv", 459e-9, 3.645409e-21),
        "sigmaemi": readCsvLL("../csv/sigmaemi_C6.csv", 501e-9, 6.020300e-21)
    }
};

earthDopants = {
    "AC46": {
        "tauT": 1.33e-4,
        "tauD": 5.33e-4,
        "wTD": 8e8,
        "wDT": 2e8,
        "sigmaabs": readCsvLL("../csv/sigmaabs_AC46.csv", 379e-9, 1e-20),
        "sigmaemi": readCsvLL("../csv/sigmaemi_AC46.csv", 614e-9, 1e-21)
    }
};

media = {
    "air": {
        "n": lambda ll: 1
    },
    "clad": {
        "n": lambda ll: 1.4,
        "att": readCsvLL("../csv/attenuationPMMA.csv", 500e-9, 0.01842*0.5)
    },
    "PMMA": {
        "n": lambda ll: math.sqrt(1/(-8.226e-15/ll**2+0.8393)+1),
        "att": readCsvLL("../csv/attenuationPMMA.csv", 500e-9, 0.01842)
    }
};

lambdaRanges = {
    "Rh6G": [410e-9, 710e-9],
    "RhB": [360e-9, 760e-9],
    "LumR": [360e-9, 750e-9],
    "LumO": [320e-9, 740e-9],
    "LumY": [360e-9, 650e-9],
    "C1": [270e-9, 680e-9],
    "C6": [320e-9, 680e-9],
    "AC46": [250e-9, 750e-9]
};
