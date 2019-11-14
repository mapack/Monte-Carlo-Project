import numpy as np

def sky(L,pos,tol):
    #Equation for the sky  
    R = L/2
    r = np.sqrt((pos[0]-R)**2+(pos[1]-R)**2+(pos[2]-R)**2)
    if np.abs(R-r) < tol:
        return 1
    else: 
        return 0

def matSci(lmbda,mat):
    #Returns the sigma, albedeo, and cos(theta) grain values based on the wavelength and materials (*no material dependence yet*)
    #Uses linear interpolation to return parameters
    data = np.loadtxt('dustdata.txt')
    lam = data[:,0]
    albedogrid = data[:,1]
    ggrid = data[:,2]
    
    falb = interpolate.interp1d(lam,albedogrid)
    albedo = falb(lmbda)
    
    fg = interpolate.interp1d(lam,ggrid)
    g = fg(lmbda)
    
    sigma = 0.

    return sigma,albedo,g

def genDirection(khat,g):
    theta = 0.0
    phi = 0.0
    
    p = np.random.rand(2)
    theta =  ((1+g)**2 - ((1-g**2)/(1-g+2*g*p[0]))**2)/(2*g)
    phi = 2*np.pi*p[1]
    k = np.zeros(2)
    k[0] = theta
    k[1] = phi
    
    return k

def posUpdate(pos,khat,step):
    theta = khat[0]
    phi = khat[1]
    npos = np.zeros(pos.shape)
    npos[0] = pos + 0 #fix later
    npos[1] = pos + 0
    npos[2] = pos + 0
    return npos

def odsSample(pos,sigma,omega,g,density,tol,Kobs,**kwargs):
    oda = 0.0 
    for key in kwargs:
        if (key == 'oda'):
            oda = kwargs[key]
    
    C = density[0][0,0,:].size
    L = density[1][-1]
    khat = Kobs
    
    if sky(pos,L,tol):
        return oda
    else:
        p = np.random.rand(1)    
        ods = -np.log(p)
        odsp = 0.0
        
        #obtain local density 
        location = np.histogramdd(pos,bins = (C,C,C),range=[(0,L),(0,L),(0,L)])
        indx = np.transpose(np.nonzero(location[0]))
        locDen = density[indx[0],indx[1],indx[2]]
        
        l = ods / (sigma * locDen)
        lcube = 0 #write a function for this 
        
        if lcube < l:
            
        if lcube > l:
            
        if lcube == l:
            odsp = sigma * locDen * l
        #odsp += sigma * locDen * step #step is the maginitude of the displacement in khat direction
        #pos += step  #needs to be vectorized to incriment in the current direction khat
        
        oda += (1/omega - 1)*odsp
        khat = genDirection(khat,g)
        
        odsSample(pos,sigma,omega,g,density,tol,khat,oda = oda)
        

def monteCarlo(density,sig,mat,lmbda,Aobs,Kobs,M,tol):
    sigma,omega,g = matSci(lmbda,mat)
    K = Kobs.size
    intensities = np.zeros(Aobs.size)        
    
    for a in range(Aobs.size):
        W = np.zeros([K,M])
        for k in range(K): 
            for m in range(M):
                pos = Aobs[a]
                oda = odsSample(pos,sigma,omega,g,density,tol,Kobs[k])
                W[k,m] = oda    
                    
        intensities[a] = 1/(K*M)*np.sum(W)           
    
    return intensities

