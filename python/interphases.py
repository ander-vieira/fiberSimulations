import math;
from math import sqrt;
from random import random;

import vector3;

# These classes represent different surfaces, used for mirrors, refractive surfaces and such
# All the classes define some common methods, which make them interchangeable and
# transparent to outside classes:
#
# intersectionPoint(pos, vel): takes a ray's position and velocity, and returns the
# distance ds such that pos+vel*ds is the first point in the surface the ray will
# reach. If no point will ever be reached (with POSITIVE ds), returns infinity.
#
# getNormalVector(pos): pos should be a position on the surface (obtained via intersectionPoint),
# and this method returns the normal vector to the surface at that point

# Describes a plane surface in 3D space
class PlaneInterphase:
    # normalVector: vector perpendicular to the plane at all points
    # origin: any point that is supposed to be in the plane
    def __init__(self, normalVector, origin=vector3.new(0, 0, 0)):
        self.normalVector = normalVector;
        self.origin = origin;

    def intersectionPoint(self, pos, vel):
        posR = pos-self.origin;

        a = vector3.dot(vel, self.normalVector);
        b = vector3.dot(posR, self.normalVector);
        if a == 0:
            ds = math.inf;
        else:
            ds = -b/a;

        if ds <= 0:
            ds = math.inf;

        return ds;

    def getNormalVector(self, pos):
        return self.normalVector;

# Describes a cylindrical surface in 3D space
class CylinderInterphase:
    # R: radius of the cylinder
    # origin: a point in the CENTER of the cylinder
    # direction: the direction of the cylinder's axis (Z axis by default)
    def __init__(self, R, origin=vector3.new(0, 0, 0), direction=vector3.new(0, 0, 1)):
        self.R2 = R**2;
        self.origin = origin;
        self.direction = direction;

    def intersectionPoint(self, pos, vel):
        posR = pos-self.origin;

        # Calculate dot products and norms using the vector's component that is
        # perpendicular to the cylinder's axis (direction).
        # These are algebraically simplified versions for better performance
        posDotDir = vector3.dot(posR, self.direction);
        velDotDir = vector3.dot(vel, self.direction);
        posDotVel = vector3.dot(posR, vel)-posDotDir*velDotDir;
        normPos2 = vector3.norm2(posR)-posDotDir**2;
        normVel2 = vector3.norm2(vel)-velDotDir**2;

        if normVel2 == 0:
            return math.inf;

        A = posDotVel/normVel2;
        B = (normPos2-self.R2)/normVel2;
        C = A**2-B;

        if C <= 0:
            return math.inf;

        if B >= 0:
            if A >= 0:
                ds = math.inf;
            else:
                ds = -sqrt(C)-A;
        else:
            ds = sqrt(C)-A;

        if ds <= 0:
            ds = math.inf;

        return ds;

    # Normal vector is defined as pointing OUTSIDE the cylinder
    def getNormalVector(self, pos):
        posR = pos-self.origin;
        result = vector3.projectTangent(posR, self.direction);

        return result/vector3.norm(result);

# Describes a parabolic surface in 3D space
class ParabolaInterphase:
    # focalPoint: distance of the focal point to the base/minimum of the parabola
    # origin: a point at the BASE of the parabola
    # directionX: the equivalent X axis in a 2D parabola (y=x^2)
    # directionY: the equivalent Y axis in a 2D parabola (y=x^2)
    # directionX and directionY should be perpendicular, otherwise bad things might happen
    def __init__(self, focalPoint, origin=vector3.new(0, 0, 0), directionX=vector3.new(1, 0, 0), directionY=vector3.new(0, 1, 0)):
        self.C = 1/4/focalPoint;
        self.origin = origin;
        self.directionX = directionX;
        self.directionY = directionY;

    def intersectionPoint(self, pos, vel):
        posR = pos - self.origin;

        posX = vector3.dot(posR, self.directionX);
        posY = vector3.dot(posR, self.directionY);
        velX = vector3.dot(vel, self.directionX);
        velY = vector3.dot(vel, self.directionY);

        if velX == 0:
            if velY == 0:
                return math.inf;
            else:
                ds = (posX**2*self.C-posY)/velY;

                if ds <= 0:
                    return math.inf;
                else:
                    return ds;

        A = (posX*velX-velY/2/self.C)/velX**2;
        B = (posX**2-posY/self.C)/velX**2;
        C = A**2-B;

        if C <= 0:
            return math.inf;

        if B >= 0:
            if A >= 0:
                ds = math.inf;
            else:
                ds = -sqrt(C)-A;
        else:
            ds = sqrt(C)-A;

        if ds <= 0:
            ds = math.inf;

        return ds;

    # Normal vector is defined as pointing UPWARDS from the parabola i.e. into
    # the cavity
    def getNormalVector(self, pos):
        posR = pos - self.origin;

        posX = vector3.projectNormal(posR, self.directionX);

        result = -2*self.C*posX+self.directionY;

        return result/vector3.norm(result);
