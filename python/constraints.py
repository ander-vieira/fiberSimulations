import vector3;

def Interval(Zmin, Zmax, reverse = False, axis = 2):
    return lambda pos: (pos[axis] >= Zmin and pos[axis] < Zmax)^reverse;

def Cylinder(R, origin = vector3.new(0, 0, 0), reverse = False):
    def result(pos):
        posR = pos-origin;

        return (posR[0]*posR[0]+posR[1]*posR[1] - R*R < 0)^reverse;

    return result;

def AlwaysTrue():
    return lambda pos: True;

def join(cond1, cond2):
    return lambda pos: (cond1(pos) and cond2(pos));
