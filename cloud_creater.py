import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations

class Cloud:
    def __init__(self, D, N, L, C, uniform=False, full_cloud=False):
        self.D = D
        self.N = N
        self.L = L
        self.C = C
        
        if uniform:
            self.point_array = generateUniformCloud(D,N,L,C)
            #self.density_array = densityFunc(self.point_array,L,C)
        if not uniform:
            self.point_array = generateCloud(D,N,L,C)
            #self.density_array = densityFunc(self.point_array,L,C)
        if full_cloud:
            self.point_array = generateFullCloud(D,N,L,C)
            #self.density_array = densityFunc(self.point_array,L,C)

def hardWallBound(arr,L):
    x = arr[:,0]
    y = arr[:,1]
    z = arr[:,2]
    
    for i in range(x.size):
        xi = x[i]
        yi = y[i]
        zi = z[i]
        
        if xi < 0:
            x[i] = xi + L
        if xi > L:
            x[i] = xi - L
        if yi < 0:
            y[i] = yi + L
        if yi > L:
            y[i] = yi - L
        if zi < 0:
            z[i] = zi + L
        if zi > L:
            z[i] = zi - L
    
    cloud = np.zeros([x.size,3])
    cloud[:,0] = x
    cloud[:,1] = y
    cloud[:,2] = z
    
    return cloud

def densityFunc(arr,L,C):
    #C is the number of cubes per axis length L
    x,y,z = np.squeeze(arr[:,0]),np.squeeze(arr[:,1]),np.squeeze(arr[:,2])
    density = np.zeros([C,C,C])
    
    boundx = np.linspace(0.0,L,C+1)
    boundy = np.linspace(0.0,L,C+1)
    boundz = np.linspace(0.0,L,C+1)
    #print(boundx)
    
    for k in range(C):
        zcondb = z > boundz[k] 
        zconda = z < boundz[k+1]
        zcond = np.logical_and(zcondb,zconda)
        for i in range(C):
            xcondb = x > boundx[i] 
            xconda = x < boundx[i+1]
            xcond = np.logical_and(xcondb,xconda)
            for j in range(C):
                ycondb = y > boundy[j] 
                yconda = y < boundy[j+1]
                ycond = np.logical_and(ycondb,yconda)
                
                #print(np.all([xcond,ycond,zcond]))
                cond = np.logical_and(np.logical_and(xcond,ycond),zcond)
                #print(cond)
                density[i,j,k] = np.sum(cond)
    
    return density   

def generateUniformCloud(D,N,L):
    return np.zeros(0)

def generateFullCloud(D,N,L):
    return np.zeros(0)

def generateCloud(D,N,L,C):
    __p1__ = np.zeros([N,3])
    __p2__ = np.zeros([N**2,3])
    __p3__ = np.zeros([N**3,3])
    __p4__ = np.zeros([N**4,3])
    delta = np.exp(np.log(N)/D)
    print("Delta: " + str(delta))
    
    for i in range(N):
        __r1__ = np.random.uniform(0,1,3)
        __p1__[i,:] = (L)*__r1__
                
    for n in range(N):
        __origin__ = np.copy(__p1__[n,:])
        #print(origin)
        for m in range(N):
            __r2__ = np.random.uniform(-1,1,3)
            __p2__[m+n*N, :] = (L/(2*delta))*__r2__ + __origin__
                
    for k in range(N**2):
        __origin__ = np.copy(__p2__[k,:])
        #print(origin)
        for q in range(N):
            __r3__ = np.random.uniform(-1,1,3)
            __p3__[q+k*N, :] = (L/(2*delta))*__r3__ + __origin__
            
    for l in range(N**3):
        __origin__ = np.copy(__p3__[l,:])
        #print(origin)
        for h in range(N):
            __r4__ = np.random.uniform(-1,1,3)
            __p4__[h+l*N, :] = (L/(2*delta))*__r4__ + __origin__
            
    
    cloud = np.concatenate((__p1__, __p2__, __p3__, __p4__), axis=0)
    cloudbn = hardWallBound(cloud,L)
    
    #density = densityFunc(cloudbn,L,C)
    density = np.histogramdd(cloudbn,bins = (C,C,C),range=[(0,L),(0,L),(0,L)])
    print("The cluster lentgh: " + str(L/(2*delta)))
    
    return cloudbn , density

D = 2.6
N = 5
L = 64
C = 64

cloudR = Cloud(D,N,L,C).point_array
cloud = cloudR[0]
print('Density Array: ' + str(cloudR[1][0]))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(cloud[:,0], cloud[:,1], cloud[:,2], c='r', marker='.')
r = [0, L]
for s, e in combinations(np.array(list(product(r, r, r))), 2):
    if np.sum(np.abs(s-e)) == r[1]-r[0]:
        ax.plot3D(*zip(s, e), color="b")

plt.show()
