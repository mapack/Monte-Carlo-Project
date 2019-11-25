import numpy as np

class Cloud:
    def __init__(self, D, N, L, C, uniform=False, full_cloud=False):
        self.D = D
        self.N = N
        self.L = L
        self.C = C
        
        if uniform:
            self.point_array = generateUniformCloud(D,N,L,C)
        if not uniform:
            self.point_array = generateCloud(D,N,L,C)
        if full_cloud:
            self.point_array = generateFullCloud(D,N,L,C)

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

def generateUniformCloud(D,N,L,C):
    dummy = np.zeros([1,3])
    density = np.zeros([C,C,C]) + 0.019
    return dummy , density

def generateFullCloud(D,N,L,C):
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
        for m in range(N):
            __r2__ = np.random.uniform(-1,1,3)
            __p2__[m+n*N, :] = (L/(2*delta))*__r2__ + __origin__
                
    for k in range(N**2):
        __origin__ = np.copy(__p2__[k,:])
        for q in range(N):
            __r3__ = np.random.uniform(-1,1,3)
            __p3__[q+k*N, :] = (L/(2*delta))*__r3__ + __origin__
            
    for l in range(N**3):
        __origin__ = np.copy(__p3__[l,:])
        for h in range(N):
            __r4__ = np.random.uniform(-1,1,3)
            __p4__[h+l*N, :] = (L/(2*delta))*__r4__ + __origin__
            
    
    cloud = np.concatenate((__p1__, __p2__, __p3__, __p4__), axis=0)
    cloudbn = hardWallBound(cloud,L)
    
    density = np.histogramdd(cloudbn,bins = (C,C,C),range=[(0,L),(0,L),(0,L)])
#    print("The cluster lentgh: " + str(L/(2*delta)))
    
    return cloudbn , density[0]
