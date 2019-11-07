import numpy as np

def sky(L,pos):
    #Equation for the sky  
    return 0 

def matSci(lmbda,mat):
    #Returns the sigma and albedeo grain values based on the wavelength and materials
    return 0

def genDirection(khat,g):
    theta = 0.0
    phi = 0.0
    
    p = np.random.rand(2)
    theta =  ((1+g)**2 - ((1-g**2)/(1-g+2*g*p[0]))**2)/(2*g)
    phi = 2*np.pi*p[1]
    k = khat * theta * phi # APPLY the new directions
    
    return k

def odsSample(pos,sigma,omega,g,density,tol,step,Kobs):
    oda = 0.0 
    C = density[0][0,0,:].size
    L = density[1][-1]
    khat = Kobs
    
    while sky(L,pos): #a func that returns FALSE if breaks the sphere
        p = np.random.rand(1)
        ods = -np.log(p)
        odsp = 0.0
        
        while np.abs(ods - odsp) <= tol:
            
            location = np.histogramdd(pos,bins = (C,C,C),range=[(0,L),(0,L),(0,L)])
            indx = np.transpose(np.nonzero(location))
            locDen = density[indx[0],indx[1],indx[2]]
            odsp += sigma * locDen * step #step is the maginitude of the displacement in khat direction
            pos += step  #needs to be vectorized to incriment in the current direction khat
        
        oda += (1/omega - 1)*odsp
        khat = genDirection(khat,g)
        
    return oda

def monteCarlo(density,sig,mat,lmbda,Aobs,Kobs,M,step,tol):
    sigma,omega,g = matSci(lmbda,mat)
    K = Kobs.size
    intensities = np.zeros(Aobs.size)        
    
    for a in range(Aobs.size):
        W = np.zeros([K,M])
        for k in range(K): 
            for m in range(M):
                pos = Aobs[a]
                oda = odsSample(pos,sigma,omega,g,density,tol,step,Kobs[k])
                W[k,m] = oda    
                    
        intensities[a] = 1/(K*M)*np.sum(W)           
    
    return intensities

