import math;
from math import sqrt, sin, cos, log, copysign;
from random import random;

import vector3;
import basics;
import params;
import constraints;

#######
# Generators
#######

# Generators are used to generate a ray's initial conditions: position, velocity,
# wavelength and power. They represent the sun / a solar simulator / a laser?

# Represents sunlight impinging on the system.
# A surface is taken, and the light is distributed evenly in it.
# Wavelength is generated using a given distribution (for example, AM 1.5)
# Each ray can be given a fraction of the incident sunlight
class SolarGenerator:
    # solarFunc: the wavelength spectrum to generate
    # ll: array of discretized wavelengths
    # direction: direction of the incident sunlight
    # vec1: one of the sides of the generator's surface
    # vec2: the other side of the generator's surface
    # origin: the origin point of the generator's surface
    # direction, vec1 and vec2 should be perpendicular
    # norm(vec1)*norm(vec2) is the surface exposed to sunlight
    def __init__(self, solarFunc, ll, direction, vec1, vec2, origin=vector3.new(0, 0, 0)):
        dLL = ll[1]-ll[0];
        self.solarDist = list(map(lambda ll: solarFunc(ll)*dLL, ll));
        self.solarConstant = sum(self.solarDist);
        self.direction = direction;
        self.vec1 = vec1;
        self.vec2 = vec2;
        self.incomingPower = self.solarConstant*vector3.norm(vec1)*vector3.norm(vec2);
        self.origin = origin;
        self.generatedPower = [0]*len(ll);

    # Generate a wavelength using the given discretized distribution
    # Returns the chosen wavelength's index k in the ll array
    def generateLambda(self):
        rand1 = random()*self.solarConstant;

        distSum = 0;
        for k in range(len(self.solarDist)):
            distSum += self.solarDist[k];
            if rand1 < distSum:
                return k;

    # Generate and return a ray with this generator's parameters:
    # Its position is uniformly distributed in the generator's surface
    # Its direction is the one given
    # Its wavelength is generated using the given distribution
    # Its power is a fraction of the incoming sunlight's power
    def generateRay(self, fraction):
        pos = self.origin+random()*self.vec1+random()*self.vec2;
        vel = self.direction;

        k = self.generateLambda();

        power = self.incomingPower*fraction;

        self.generatedPower[k] += power;

        return basics.Ray(pos, vel, k, power);

#######
# Optical components
#######

# Optical components have an effect on the light rays that pass through them.
# To be interchangeable and transparent to outside classes, they define several
# common methods:
#
# intersectionPoint: takes a ray, and returns the distance ds the ray will go
# before being affected by the optical component.
# For surface components like mirrors, this hooks into the component's interphase object.
# For medium components like attenuation, it generates an absorption/emission/whatever
# length as a random distributed variable.
#
# process: takes a ray and the distance ds from intersectionPoint.
# After the ray is moved the distance ds given by intersectionPoint, the
# component has its effect on the ray (reflection, refraction, absorption, etc.)

# A mirror (reflects all incident light)
class Mirror:
    # Interphase: an object of any interphase class, defines the component's geometry
    # cint: a constraint for the component to work
    def __init__(self, interphase, cint = constraints.AlwaysTrue()):
        self.interphase = interphase;
        self.constraint = cint;

    def setConstraint(self, func):
        self.constraint = func;

    def reset(self):
        pass;

    def intersectionPoint(self, ray):
        return self.interphase.intersectionPoint(ray.pos, ray.vel);

    # The ray's position should be on the mirror's surface
    # Reflects the ray fully off the surface if the constraint is met
    def process(self, ray, ds):
        if self.constraint(ray.pos):
            normalVector = self.interphase.getNormalVector(ray.pos);

            # Reflected direction
            ray.vel = ray.vel-2*vector3.projectNormal(ray.vel, normalVector);

# A refractive surface, such as between the air and the fiber's core
# Uses Snell's law to calculate refraction, partial reflection or total internal reflection
# depending on the incidence angle, refraction indices etc.
class Refractor:
    # Interphase: an object of any interphase class, defines the component's geometry
    # mediumPlus: the name of the outside medium as defined in params.py
    # mediumMinus: the name of the inside medium as defined in params.py
    # ll: array of discretized wavelengths
    # cint: a constraint for the component to work
    # The outside medium is the one the normal vector points at.
    def __init__(self, interphase, mediumPlus, mediumMinus, ll, cint = constraints.AlwaysTrue()):
        self.interphase = interphase;
        self.Nplus = list(map(params.media[mediumPlus]["n"], ll));
        self.Nminus = list(map(params.media[mediumMinus]["n"], ll));
        self.constraint = cint;

    def setConstraint(self, func):
        self.constraint = func;

    def reset(self):
        pass;

    def intersectionPoint(self, ray):
        return self.interphase.intersectionPoint(ray.pos, ray.vel);

    # The ray's position should be on the refractive surface
    # Applies Snell's law, calculates refraction, partial reflection
    # or TIR as needed
    def process(self, ray, ds):
        if self.constraint(ray.pos):
            normalVector = self.interphase.getNormalVector(ray.pos);

            # Get the reflected direction (for total or partial reflection)
            reflectedVel = ray.vel-2*vector3.projectNormal(ray.vel, normalVector);

            # Get the incidence angle
            cosTheta = vector3.dot(ray.vel, normalVector);

            # Choose the incident and refracted media based on the angle
            if cosTheta > 0:
                N1 = self.Nminus[ray.k];
                N2 = self.Nplus[ray.k];
            else:
                N1 = self.Nplus[ray.k];
                N2 = self.Nminus[ray.k];

            if N2 < N1 and cosTheta < sqrt(1-pow(N2/N1, 2)):
                # Total internal reflection

                ray.vel = reflectedVel;
            else:
                # Refraction or partial reflection

                # Get the refraction angle for calculating Fresnel reflection
                newCosTheta = copysign(sqrt(1-(1-cosTheta*cosTheta)*(N1/N2)**2), cosTheta);

                # Get the Fresnel coefficient for partial reflection
                fresnelR = basics.fresnelR(cosTheta, newCosTheta, N1, N2);

                if random() < fresnelR:
                    # Partial reflection (random chance)

                    ray.vel = reflectedVel;
                else:
                    # Refraction (if no partial reflection)

                    # Calculate refracted direction
                    tangentVel = ray.vel-cosTheta*normalVector;
                    newTangentVel = tangentVel*N1/N2;
                    newNormalVel = newCosTheta*normalVector;

                    ray.vel = newTangentVel+newNormalVel;

# A component to collect rays and add the power that reaches them
# Represents a solar cell or other kind of detector
# The output power is saved for each wavelength separately.
class Collector:
    # interphase: an object of any interphase class, defines the component's geometry
    # numLL: number of discretized wavelength values (length of ll)
    # cint: a constraint for the component to work
    def __init__(self, interphase, numLL, cint = constraints.AlwaysTrue()):
        self.interphase = interphase;
        self.finalPower = [0]*numLL;
        self.constraint = cint;
        self.counter = 0;

    def setConstraint(self, func):
        self.constraint = func;

    def reset(self):
        for k in range(len(self.finalPower)):
            self.finalPower[k] = 0;

    def getPower(self, k):
        return self.finalPower[k];

    def getTotalPower(self):
        return sum(self.finalPower);

    def intersectionPoint(self, ray):
        return self.interphase.intersectionPoint(ray.pos, ray.vel);

    # The collector only detects rays going one way (defined by the normal vector)
    def process(self, ray, ds):
        if self.constraint(ray.pos):
            normalVector = self.interphase.getNormalVector(ray.pos);

            if vector3.dot(ray.vel, normalVector) > 0:
                self.finalPower[ray.k] += ray.power;
                self.counter += 1;
                ray.loopOn = False;

# Represents attenuation caused by a medium such as PMMA
class Attenuation:
    # medium: the name of the medium that causes the attenuation (ex. "PMMA")
    # Defined in params.py
    # ll: array of discretized wavelengths
    # cint: a constraint for the component to work
    def __init__(self, medium, ll, cint=constraints.AlwaysTrue()):
        self.alpha = list(map(params.media[medium]["att"], ll)); # m^-1
        self.constraint = cint;

    def setConstraint(self, func):
        self.constraint = func;

    def reset(self):
        pass;

    # The point of absorption follows an exponential distribution
    def intersectionPoint(self, ray):
        if self.constraint(ray.pos):
            if self.alpha[ray.k] <= 0:
                return math.inf;
            else:
                return -log(random())/self.alpha[ray.k];
        else:
            return math.inf;

    # Absorb the ray: stop simulating it
    def process(self, ray, ds):
        if self.constraint(ray.pos):
            ray.loopOn = False;

# Represents absorption and emission of the light by a dye dopant
# Stimulated emission (and N2 effects in general) isn't considered yet
class DyeDopant:
    # dopant: the name of the dye dopant (ex. "Rh6G"), defined in params.py
    # N: the concentration of the dopant in molecules/m^3
    # ll: array of discretized wavelengths
    # cint: a constraint for the component to work
    def __init__(self, dopant, N, ll, cint=constraints.AlwaysTrue()):
        self.ll = ll;
        sigmaabsFun = params.dyeDopants[dopant]["sigmaabs"];
        sigmaemiFun = params.dyeDopants[dopant]["sigmaemi"];
        self.sigmaabs = list(map(sigmaabsFun, ll)); # m^2
        self.sigmaemi = list(map(sigmaemiFun, ll)); # m^2
        self.alpha = list(map(lambda sigma: sigma*N, self.sigmaabs));
        self.sumEmi = sum(self.sigmaemi);
        tauRad = params.dyeDopants[dopant]["tauRad"];
        tauNR = params.dyeDopants[dopant]["tauNR"];
        self.quantumYield = tauNR/(tauRad+tauNR);
        self.N = N;
        self.constraint = cint;

    def setConstraint(self, func):
        self.constraint = func;

    def reset(self):
        pass;

    # The point of absorption follows an exponential distribution
    def intersectionPoint(self, ray):
        if self.constraint(ray.pos):
            if self.alpha[ray.k] <= 0:
                return math.inf;
            else:
                return -log(random())/self.alpha[ray.k];
        else:
            return math.inf;

    # Generate a new wavelength using the emission spectrum as a distribution
    def generateLambda(self):
        rand1 = random()*self.sumEmi;

        sigmaSum = 0;
        for k in range(len(self.sigmaemi)):
            sigmaSum += self.sigmaemi[k];
            if rand1 < sigmaSum:
                return k;

    # Absorb the ray and spontaneously emit it
    # There's no time variable, so emission happens immediately after
    # Stimulated emission etc. aren't included yet.
    def process(self, ray, ds):
        if self.constraint(ray.pos):
            # Generate random direction after spontaneous emission
            rand1 = 2*math.pi*random();
            rand2 = 2*random()-1;
            ray.vel = vector3.new(sqrt(1-pow(rand2,2))*cos(rand1), sqrt(1-pow(rand2,2))*sin(rand1), rand2);

            # Generate new wavelength using emission spectrum
            oldK = ray.k;
            ray.k = self.generateLambda();

            # Change photon power (Stokes shift)
            ray.power = ray.power*self.ll[oldK]/self.ll[ray.k];

            # Quantum yield: some spontaneous emissions are lost due to
            # nonradiative decay
            if random() > self.quantumYield:
                ray.loopOn = False;

#######
# Miscellaneous components
#######

# These components don't fit in with generator or optical components
# They don't follow any particular model

# Combine several Collector components and add their outputs together
class CollectorGroup:
    def __init__(self, cllList):
        self.cllList = cllList;

    def getPower(self, k):
        result = 0;
        for cll in self.cllList:
            result += cll.getPower(k);

        return result;

    def getTotalPower(self):
        result = 0;
        for cll in self.cllList:
            result += cll.getTotalPower();

        return result;

    def reset(self):
        for cll in self.cllList:
            cll.reset();
