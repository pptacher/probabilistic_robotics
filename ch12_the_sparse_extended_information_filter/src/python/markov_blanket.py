from numpy import *

from pdb import set_trace as bp

def markov_blanket(n,m0,G):
    s = G.shape[0]
    mb1 = m0
    mb2 = G[n,:].nonzero()[1]

    mb12 = intersect_int(mb1,mb2,s)
    if mb12.size == 0:
        mb = bfs(n,G)
        mb2 = union_int(mb2,mb,s)

    mb2 = hstack(([0,n],mb2))
    mb = union_int(mb1,mb2,s)
    return mb

def union_int(m1,m2,s):
    x = zeros(s)
    y = zeros(s)
    x[m1] = 1
    y[m2] = 1
    return nonzero(x+y)[0]

def intersect_int(m1,m2,s):
    x = zeros(s)
    y = zeros(s)
    x[m1] = 1
    y[m2] = 1
    return nonzero(x*y)[0]

def bfs(n,G):
    s = G.shape[0]
    d = -ones(s,dtype=int)
    d[0] = 0
    q = zeros(s,dtype=int)
    i=0
    j=1
    q[0]=0
    found=False
    p=zeros(0)
    curr=0
    while (not found) and i<j:
        curr = q[i]
        i += 1
        next = nonzero(logical_and(G[curr,:].astype(int).toarray().ravel(),d==-1))[0]
        s1 = next.size
        q[j:j+s1]=next
        j += s1
        d[next]=curr
        found=d[n]==curr

    if found:
        p = curr
        prev = d[curr]
        while prev!=0:
            p = hstack((p,prev))
            prev = d[prev]

    return p
