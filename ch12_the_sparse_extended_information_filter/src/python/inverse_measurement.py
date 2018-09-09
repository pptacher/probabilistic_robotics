from numpy import *

def inverse_measurement(m,z):
    c = cos(m[2]+z[1,:])
    s = sin(m[2]+z[1,:])

    return m[0:2].reshape((2,1)) + vstack((z[0,:]*s,-z[0,:]*c))
