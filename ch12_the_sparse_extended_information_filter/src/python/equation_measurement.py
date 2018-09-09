from numpy import *

def equation_measurement(m,ind):
    N = ind.size
    ind1 = 1+2*ind
    x = m[ind1]
    y = m[ind1+1]

    delta = vstack((x,y))-m[0:2].reshape((2,1))
    q = sqrt(sum(delta**2,axis=0))
    return vstack((q,arctan2(delta[1,:],delta[0,:])-m[2]+pi/2))
