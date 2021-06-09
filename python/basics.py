class Ray:
    def __init__(self, pos, vel, k, power):
        self.loopOn = True;
        self.pos = pos;
        self.vel = vel;
        self.k = k;
        self.power = power;

    def move(self, ds):
        self.pos += self.vel*ds;

# Return the Fresnel reflection coefficient for given incidence angle and
# refraction indices.
def fresnelR(cos1, cos2, N1, N2):
    return (pow((N1*cos1-N2*cos2)/(N1*cos1+N2*cos2), 2)+pow((N2*cos1-N1*cos2)/(N2*cos1+N1*cos2), 2))/2;
