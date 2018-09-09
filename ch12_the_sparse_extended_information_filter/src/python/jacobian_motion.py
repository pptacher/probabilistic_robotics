from numpy import *

def jacobian_motion(theta,v,a,dt):
    H = 0.74
    L = 2.83
    a1 = 0.95 + L
    b1 = 0.5

    return array([[1,0,-dt*(v*sin(theta)+v/L*tan(a)*(a1*cos(theta)-b1*sin(theta)))], \
                           [0,1,dt*(v*cos(theta)-v/L*tan(a)*(a1*sin(theta)+b1*cos(theta)))], \
                           [0,0,1]])
