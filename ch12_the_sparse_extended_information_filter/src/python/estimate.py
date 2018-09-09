from numpy import *
from scipy.linalg import solve

from pdb import set_trace as bp

def estimate(m0,xi,omega,m):
    n = (xi.size-3)//2

    if n < 100:

       m1 = solve(omega,xi)

    else:
       m1=m
       for i in range(0,m0.size):
           ind = 2*m0[i]+1

           m1[ind:ind+2] = solve(omega[ind:ind+2,ind:ind+2],xi[ind:ind+2]-omega[ind:ind+2,:].dot(m1)+omega[ind:ind+2,ind:ind+2].dot(m1[ind:ind+2]))
       m1[0:3] = solve(omega[0:3,0:3],xi[0:3]-omega[0:3,:].dot(m1)+omega[0:3,0:3].dot(m1[0:3]))
    return m1
