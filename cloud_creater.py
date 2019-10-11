import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def cloud_creator(D,N,L):
    p1 = np.zeros([N,3])
    p2 = np.zeros([N**2,3])
    p3 = np.zeros([N**3,3])
    p4 = np.zeros([N**4,3])
    delta = np.exp(np.log(N)/D)
    print("Delta: " + str(delta))
    
    for i in range(N):
        r1 = np.random.uniform(-1,1,3)
        p1[i,:] = (L/2)*r1
                
    for n in range(N):
        origin = np.copy(p1[n,:])
        #print(origin)
        for m in range(N):
            r2 = np.random.uniform(-1,1,3)
            p2[m+n*N, :] = (L/(2*delta))*r2 + origin
                
    for k in range(N**2):
        origin = np.copy(p2[k,:])
        #print(origin)
        for q in range(N):
            r3 = np.random.uniform(-1,1,3)
            p3[q+k*N, :] = (L/(2*delta))*r3 + origin
            
    for l in range(N**3):
        origin = np.copy(p3[l,:])
        #print(origin)
        for h in range(N):
            r4 = np.random.uniform(-1,1,3)
            p4[h+l*N, :] = (L/(2*delta))*r4 + origin
            
            
    print("The cluster lentgh: " + str(L/(2*delta)))
    return (p1,p2,p3,p4)

#cloud_creater(2.6,5,1)

def main():
    D = 2.6
    N = 2
    L = 100
    data = cloud_creator(D,N,L)
    p1 = data[0]
    p2 = data[1]
    p3 = data[2]
    p4 = data[3]
    
    #print(p1)
    #print(p2)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x1 = p1[:,0]
    y1 = p1[:,1]
    z1 = p1[:,2]
    
    #print(x1)
    
    x2 = p2[:,0]
    y2 = p2[:,1]
    z2 = p2[:,2]
    
    x3 = p3[:,0]
    y3 = p3[:,1]
    z3 = p3[:,2]
    
    x4 = p4[:,0]
    y4 = p4[:,1]
    z4 = p4[:,2]
    
    ax.scatter(x1, y1, z1, c='r', marker='o')
    ax.scatter(x2, y2, z2, c='b', marker='o')
    ax.scatter(x3, y3, z3, c='m', marker='o')
    ax.scatter(x4, y4, z4, c='g', marker='o')


    plt.show()
    
main()