from numpy import *
from scipy.sparse import csr_matrix,csc_matrix,lil_matrix
from scipy.linalg import inv,solve
import numpy as np

from jacobian_motion import jacobian_motion
from equation_motion import equation_motion

from pdb import set_trace as bp

def motion(v,a,dt,m0,xi,omega,m,G):
    n = xi.size
    R = diag([8e-4,8e-4,1e-6])
    J = jacobian_motion(m[2],v,a,dt)

    psi = inv(J)-eye(3)
    ind = 2*sort(m0) + 1
    ind = vstack((ind,ind+1)).flatten(order='F')
    ind  = hstack((r_[0:3],ind)).astype(int)
    omega_s = omega[ix_(ind,r_[0:3])]
    psi1 = omega_s@psi

    lmbd = zeros([n,n])
    lmbd[ix_(ind,r_[0:3])] = psi1;lmbd[ix_(r_[0:3],ind)]=lmbd[ix_(r_[0:3],ind)]+psi1.T
    lmbd[0:3,0:3] =  lmbd[0:3,0:3] + psi.T@omega[0:3,0:3]@psi
    phi = omega + lmbd

    phiX = phi[ix_(r_[0:3],ind)]
    k=zeros([n,n])
    k[ix_(ind,ind)] = phiX.T@inv(inv(R)+phi[0:3,0:3])@phiX
    omega1 = phi-k

    dx = equation_motion(m[2],v,a,dt)
    omega_s = omega1[ix_(ind,r_[0:3])]

    xi1 = xi
    xi1[ind] = xi1[ind] + omega_s@dx
    xi1[ind] = xi1[ind] + (lmbd-k)[ind,:]@m
    m1 = m
    m1[0:3] = m1[0:3] + dx

    if m0.size>0:
        s,t = meshgrid(m0,m0)
        s = s.flatten()
        t = t.flatten()
        edges = sort(np.vstack((s,t)),axis=0)
        edges = unique(edges,axis=1)   #needs version >=1.13
        #edges = np.vstack({tuple(row) for row in edges.T}).T
        edges = edges[:,edges[0,:]-edges[1,:]!=0]
        if edges.shape[1]>0:
            sg = G.shape[0]

            nl = csr_matrix((ones(edges.shape[1]),(edges[0,:],edges[1,:])),shape=(sg,sg),dtype=bool)
            G = G + nl + nl.T

    return xi1,omega1,m1,G
