from numpy import *
import numpy as np
from scipy.linalg import inv
from scipy.sparse import csr_matrix,vstack,hstack

from jacobian_measurement import jacobian_measurement
from equation_measurement import equation_measurement
from inverse_measurement import inverse_measurement

from pdb import set_trace as bp

def measurement(z,c,xi,omega,m,G):
    n = (xi.size-3)//2
    Q = diag([5.0,0.02])
    c = c.astype(int)
    new = max(c)-n

    if new>0:
        nsize = xi.size+2*new
        omega1 = zeros([nsize,nsize])
        omega1[0:xi.size,0:xi.size] = omega
        xi1 = zeros(nsize)
        xi1[0:xi.size] = xi
        m1 = zeros(nsize)
        m1[0:xi.size] = m
        m1[xi.size:] = inverse_measurement(m[0:3],z[:,c>n]).T.ravel()
        s = G.shape[0]
        G.resize((s+new,s+new))
    else:
        omega1=omega
        xi1=xi
        m1=m

    J = jacobian_measurement(m1,c)
    zhat = equation_measurement(m1,c)

    for i in range(0,c.size):
           j = c[i]
           h = J[:,:,i]
           dz = z[:,i]-zhat[:,i]

           dz[1] = measure(dz[1])
           dz = dz.reshape((2,1))
           ind = hstack((r_[0:3],[2*j+1,2*j+2])).toarray().ravel()

           Q1 = inv(Q)
           xi1[ind] = xi1[ind]+ h.T@Q1@(dz+h@m1[ind].reshape((5,1))).flatten()
           omega1[ix_(ind,ind)] = omega1[ix_(ind,ind)] + h.T@Q1@h

    c = unique(c)
    row = zeros(c.size)
    col = c
    data = ones(c.size)
    nl = csr_matrix((data,(row,col)),shape=(G.shape[0],G.shape[0]),dtype=bool)
    G = G + nl + nl.T

    return xi1,omega1,m1,G


def measure(theta):
                   tmp = theta%(2*pi)
                   return tmp-minimum(floor(tmp/pi),1)*2*pi
