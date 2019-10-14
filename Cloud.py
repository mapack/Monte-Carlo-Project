import numpy as np

class Cloud:
    def __init__(self, D, N, L, uniform=False, full_cloud=False):
        self.D = D
        self.N = N
        self.L = L
        
        if uniform:
            self.point_array = generateUniformCloud(D,N,L)
            self.density_array = np.zeros(0)
        if not uniform:
            self.point_array = generateCloud(D,N,L)
            self.density_array = np.zeros(0)
        if full_cloud:
            self.point_array = generateFullCloud(D,N,L)
            self.density_array = np.zeros(0)

   

def generateUniformCloud(D,N,L):
    return np.zeros(0)

def generateFullCloud(D,N,L):
    return np.zeros(0)

def generateCloud(D,N,L):
    __p1__ = np.zeros([N,3])
    __p2__ = np.zeros([N**2,3])
    __p3__ = np.zeros([N**3,3])
    __p4__ = np.zeros([N**4,3])
    delta = np.exp(np.log(N)/D)
    print("Delta: " + str(delta))
    
    for i in range(N):
        __r1__ = np.random.uniform(-1,1,3)
        __p1__[i,:] = (L/2)*__r1__
                
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
            
    print("The cluster lentgh: " + str(L/(2*delta)))
    return np.concatenate((__p1__, __p2__, __p3__, __p4__), axis=0)
