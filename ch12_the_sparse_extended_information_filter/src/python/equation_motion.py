from numpy import *

def equation_motion(theta,v,a,dt):
    H = 0.74
    L = 2.83
    a1 = 0.95 + L
    b1 = 0.5

    v = v/(1-tan(a)*H/L)

    return array([dt*(v*cos(theta)-v/L*tan(a)*(a1*sin(theta)+b1*cos(theta))),\
                  dt*(v*sin(theta)+v/L*tan(a)*(a1*cos(theta)-b1*sin(theta))),\
                  dt*v/L*tan(a)])
