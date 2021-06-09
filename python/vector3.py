import numpy;
from math import sqrt;

# Create a new vector from Cartesian coordinates
def new(x, y, z):
    return numpy.asarray([x, y, z]);

# Get the squared norm of a vector (avoid square root for performance)
def norm2(vec):
    return vec[0]**2+vec[1]**2+vec[2]**2;

# Get the norm of a vector
def norm(vec):
    return sqrt(norm2(vec));

# Get dot product / scalar product of two vectors
def dot(vec1, vec2):
    return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];

# Project a vector onto a normal vector
def projectNormal(vec1, normalVector):
    return dot(vec1, normalVector)*normalVector;

# Project a vector onto the surface perpendicular to a normal vector
def projectTangent(vec1, normalVector):
    return vec1-projectNormal(vec1, normalVector);
