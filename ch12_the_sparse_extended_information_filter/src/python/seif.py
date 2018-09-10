from numpy import *
import scipy.linalg as la
from scipy.sparse import csr_matrix, save_npz
import scipy.io as sio
import os
import time as timechr

import matplotlib.pyplot as plt
from motion import motion
from estimate import estimate
from correspondence import correspondence
from measurement import measurement
from sparsification import sparsification

from pdb import set_trace as bp

def seif():
    filename1 = 'z.mat'
    filename2 = 'aa3_dr.mat'
    filename3 = 'aa3_lsr2.mat'
    filename4 = 'c_5000.mat'
    include_dir = './data/'

    z_contents = sio.loadmat(include_dir +filename1)
    z = z_contents['z'];

    z_contents = sio.loadmat(include_dir +filename2)
    speed = z_contents['speed'].ravel();
    steering = z_contents['steering'].ravel();
    time = z_contents['time'].ravel();

    z_contents = sio.loadmat(include_dir +filename3)
    timeLsr = z_contents['TLsr'].ravel();
    L = size(timeLsr,0);

    z_contents = sio.loadmat(include_dir +filename4)
    corresp = z_contents['corresp']
    Lc = size(corresp,0);

    del z_contents

    dt = 25e-3
    G = csr_matrix((1,1),dtype=bool)
    Q = diag([5.0,0.02])
    m = zeros(3)
    xi = zeros(3)
    omega = 10e4*eye(3)
    m0=zeros(0)

    #max active landmarks
    N = 20
    plt.ion()

    fig,ax = plt.subplots(1,1)
    ax.set_aspect('equal')
    ax.set_xlim(-100, 300)
    ax.set_ylim(-50, 350)
    ax.hold(True)
    plt.show(False)
    plt.draw()
    #background = fig.canvas.copy_from_bbox(ax.bbox)
    line1 = ax.plot(0,0,'b-')[0]
    line2 = ax.plot(-1000,1000,'ro',markersize=2)[0]

    poses = zeros([3,5000])
    stindex=0
    j=searchsorted(timeLsr,time[stindex])
    timechr.sleep(3)
    t1 = timechr.time()

    for i in range(stindex, time.shape[0]):
    #for i in range(stindex, 10000):
        t3= timechr.time()

        if i>stindex and i%5000==0:
            save(include_dir +'poses_'+str(i),poses)
            save(include_dir +'landmarks_'+str(i),m)
            save(include_dir +'xi_'+str(i),xi)
            save(include_dir +'omega_'+str(i),omega)
            save(include_dir +'m0_'+str(i),m0)
            save_npz(include_dir +'G_', G)

        xi,omega,m,G = motion(speed[i],steering[i],dt,m0,xi,omega,m,G)

        m = estimate(m0,xi,omega,m)

        while j<L and timeLsr[j]<time[i]+dt*1000:
           
           z1=z[0:2,nonzero(z[0,:,j]),j].transpose((0,2,1))[:,:,0]
           if z1.size >0:
               co = correspondence(z1,m0,omega,m,G)
               xi,omega,m,G = measurement(z1,co,xi,omega,m,G)
           j+=1

           if co.size>0:
               n=(xi.size-3)//2
               m0a = setdiff1d(m0,co)
               #returns float when m0 empty
               tq = hstack((m0a,unique(co))).astype(int)
               tq = tq[max(0,tq.size-N):]
               m1 = setdiff1d(m0,tq)
               m1 = union1d(m1,setdiff1d(co,tq))
               m0 = tq
               if m1.size>0:
                   xi,omega,G = sparsification(m0,m1,xi,omega,m,G)

        print("\x1b[2J\x1b[H")
        t2=timechr.time()
        print('iter: '+str(i))
        print('iter time: {0:.5f}'.format(t2-t3))
        print('avg: {0:.5f}'.format((t2-t1)/(i+1)))
        poses[:,i%5000]=m[0:3]

        if i%5000==4999:
            line2.set_data(m[3::2],m[4::2])
            line1.set_xdata(hstack((line1.get_xdata(),poses[0,:])))
            line1.set_ydata(hstack((line1.get_ydata(),poses[1,:])))
            #fig.canvas.restore_region(background)
            #ax.draw_artist(line1)
            #ax.draw_artist(line2)
            fig.canvas.blit(ax.bbox)
            fig.canvas.draw()
            #plt.show()
            #fig.canvas.flush_events()


if __name__ == "__main__":
    # execute only if run as a script
    import cProfile
    pr = cProfile.Profile()
    pr.enable()
    seif()
    pr.disable()
    pr.print_stats(sort='cumtime')
