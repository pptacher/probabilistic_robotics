from numpy import *
import numpy as np
import scipy.linalg as la
from scipy.stats import chi2
from scipy.sparse import csr_matrix
from equation_measurement import equation_measurement
from jacobian_measurement import jacobian_measurement
from markov_blanket import markov_blanket

from pdb import set_trace as bp

#@profile
def correspondence(z,m0,omega,m,G):

    k = z.shape[1]
    if k==0:
        return zeros([1,0])
    n = (m.size-3)//2
    Q = diag([5.0,0.02])

    co=zeros(k,dtype=int)

    if n>0:

        t = m[3:]
        t1 = t[0::2]
        t2 = t[1::2]
        lst  = nonzero(logical_and(logical_and(abs(t1-m[0])<= 75\
                       , abs(t2-m[1])<= 75)\
                       , (t1-m[0])*cos(m[2])+(t2-m[1])*sin(m[2]) >= -5 ))[0]
        a  = lst.size
        aa = zeros([2*a,2])

        J = jacobian_measurement(m,lst+1)
        zhat = equation_measurement(m,lst+1).reshape((2,1,a))

        for i in range(0,a):
                ind = markov_blanket(lst[i]+1,m0,G)
                ind = ind[1:]-1
                c1 = 3+2*ind
                ind = vstack((c1,c1+1)).flatten(order='F')
                ind1 = nonzero(ind==3+2*lst[i])[0][0]
                f = fxmd(ind1,3+ind.size)
                ind2 = np.hstack((r_[0:3],ind))
                ind3 = np.hstack((r_[0:3],array([3+ind1,4+ind1])))

                s = la.solve(omega[ix_(ind2,ind2)],f.T@J[:,:,i].T)
                d,u = la.eigh(J[:,:,i]@s[ind3,:]+Q)
                u = diag(sqrt(1/d)).dot(u.T)
                aa[2*i:2*i+2,:] = u

        f =  np.transpose(zhat-z.reshape((2,k,1)),(2,0,1)).reshape((2*a,k))
        f[1::2,:] = measure(f[1::2,:])

        d1 = zeros([2*a,k]);
        ia = arange(a)[in1d(lst, m0-1)]
        cc = lst[ia]

        for j in range(0,cc.size):
                ind=ia[j]
                d1[2*ind:2*ind+2,:] = aa[2*ind:2*ind+2,:].dot(f[2*ind:2*ind+2,:])

        d2 = d1[vstack((2*ia,2*ia+1)).flatten(order='F'),:]
        d2 = d2**2
        m2 = d2[1::2,:]
        d3 = d2[0::2,:] + m2

        m1 = amin(d3,axis=0)
        co = argmin(d3,axis=0)
        new = nonzero(m1>chi2.ppf(0.50,2))[0]
        new1 = nonzero(m1[new]<chi2.ppf(0.95,2))[0]
        ia = ia[setdiff1d(r_[0:cc.size],co[new[new1]])]

        co = cc[co]+1

        if (new.size >0 and cc.size<a) or cc.size==0:
                if cc.size==0:
                        new = r_[0:k]
                ib = setdiff1d(r_[0:a],ia)
                lst = lst[ib]

                for j in range(0,lst.size):
                        ind = ib[j]

                        d1[ix_([2*ind,2*ind+1],new)] = aa[2*ind:2*ind+2,:].dot(f[ix_([2*ind,2*ind+1],new)])
                if ib.size>0:
                        d2 = d1[ix_(vstack((2*ib,2*ib+1)).flatten(order='F'),new)]
                        d2 = d2**2
                        m2 = d2[1::2,:]
                        d3 = d2[0::2,:] + m2
                        m1 = amin(d3,axis=0)
                        co[new] = argmin(d3,axis=0)

                        co[new] = lst[co[new]]+1
                        new = new[m1>chi2.ppf(0.95,2)]
                        new1 = zeros(0)
    else:
        new = r_[0:k]
        new1 = zeros(0)

    if new1.size>0:
            new = new[setdiff1d(r_[0:new.size],new1)]
    if new.size>0:
            co[new] = n+ r_[1:new.size+1]

    return co

def fxmd(p,n):
        if p>=0:
                res=zeros([5,n])
                res[0:3,0:3] = eye(3)
                res[3,3+p] = 1
                res[4,4+p] = 1
                return res

def measure(theta):
        tmp = theta%(2*pi)
        return tmp-minimum(floor(tmp/pi),1)*2*pi
