import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def cloud_creator(D,N,L,R):
    NR = 0
    for i in range(R):
        NR += N**(R-i)
    p = np.zeros([NR,3])
    delta = np.exp(np.log(N)/D)
    print("Delta: " + str(delta))
    
    for r in range(R):
        
        if r == 0:
            for i in range(N):
                rand = np.random.uniform(-1,1,3)
                p[i,:] = (L/2) * rand
                
        else:
            for m in range(N):
                origin = np.copy(p[m + r*,:])
                for n in range(N):
                    rand = np.random.uniform(-1,1,3)
                    p[(m+*N**r + n,:] = (L/(2*delta))*rand + origin
    
    
    #for i in range(N):
    #    r1 = np.random.uniform(-1,1,3)
    #    for j in range(3):
    #        p1[i,j] = (L/2)*r1[j]
                
    #for n in range(N):
    #    origin = np.copy(p1[n,:])
    #    #print(origin)
    #    for m in range(N):
    #        r2 = np.random.uniform(-1,1,3)
    #        for p in range(3):
    #            p2[m+n*N, p] = (L/(2*delta))*r2[p] + origin[p]
                
    
                
    print("The cluster lentgh: " + str(L/(2*delta)))
    return p

def main():
    D = 2.6
    N = 2
    L = 100
    R = 4
    
    p = cloud_creator(D,N,L,R)
    
    print(p)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x1 = p[:,0]
    y1 = p[:,1]
    z1 = p[:,2]
    
    #print(x1)
    
    ax.scatter(x1, y1, z1, c='r', marker='o')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()
    
main()