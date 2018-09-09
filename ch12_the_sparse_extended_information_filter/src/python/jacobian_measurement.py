from numpy import *

from pdb import set_trace as bp

def jacobian_measurement(m,ind):
    N = ind.size
    ind1 = 1+2*ind

    x = m[ind1]
    y = m[ind1+1]

    delta = vstack((x,y))-m[0:2].reshape((2,1))
    q = sum(delta**2,axis=0)

    J = array([[1,0],[0,-1],[0,1],[1,0]]).dot(delta)
    #ugly
    J = vstack((vstack((vstack((-J,zeros([1,N]))),-ones([1,N]))),J))
    J[0::2,:] = J[0::2,:] *1/sqrt(q)
    J[array([1,3,7,9]),:] = J[array([1,3,7,9]),:] *1/q

    return J.reshape((2,5,N), order='F')
