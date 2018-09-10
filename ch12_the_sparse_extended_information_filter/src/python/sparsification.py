from numpy import *
from scipy.sparse import csr_matrix
from scipy.linalg import solve

from pdb import set_trace as bp

def sparsification(m0,m1,xi,omega,m,G):
    sg = G.shape[0]
    ol = csr_matrix((ones(m1.size),(zeros(m1.size),m1)),shape=(sg,sg),dtype=bool)
    G = G - ol - ol.T
    G = G.astype(bool)
    s,t = meshgrid(union1d(m0,m1),m1)
    s = s.flatten()
    t = t.flatten()
    edges = sort(vstack((s,t)),axis=0)
    edges = unique(edges,axis=1)   #needs version >=1.13
    #edges = vstack({tuple(row) for row in edges.T}).T
    edges = edges[:,edges[0,:]-edges[1,:]!=0]
    nl = csr_matrix((ones(edges.shape[1]),(edges[0,:],edges[1,:])),shape=(sg,sg),dtype=bool)
    G = G + nl + nl.T

    m0 = sort(m0)
    m1 = sort(m1)
    m0 = 2*m0 + 1
    m0 = vstack((m0,m0+1)).flatten(order='F')
    m1 = 2*m1 + 1
    m1 = vstack((m1,m1+1)).flatten(order='F')

    omega2 = omega[ix_(hstack((hstack((r_[0:3],m0)),m1)),m1)]
    l1 = omega2.dot(solve(omega[ix_(m1,m1)],omega2.T))

    omega3 = omega[ix_(hstack((hstack((r_[0:3],m0)),m1)),hstack((r_[0:3],m1)))]
    l2 = omega3.dot(solve(omega[ix_(hstack((r_[0:3],m1)),hstack((r_[0:3],m1)))],omega3.T))

    omega4 = omega[:,r_[0:3]]
    l3 = omega4.dot(solve(omega[ix_(r_[0:3],r_[0:3])],omega4.T))

    omega1 = copy(omega)
    omega1[ix_(hstack((hstack((r_[0:3],m0)),m1)),hstack((hstack((r_[0:3],m0)),m1)))] -= l1
    omega1[ix_(hstack((hstack((r_[0:3],m0)),m1)),hstack((hstack((r_[0:3],m0)),m1)))] += l2
    omega1 -= l3

    omega1[ix_(r_[0:3],m1)] = 0
    omega1[ix_(m1,r_[0:3])] = 0

    xi1 = xi + (omega1-omega)@m

    return xi1,omega1,G
