import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import time

def sky(pos,L,tol):
    #Equation for the sky  
    R = float(L)/2
    r = np.sqrt((pos[0]-R)**2+(pos[1]-R)**2+(pos[2]-R)**2)
    if np.abs(R-r) < tol:
        return 0
    elif r > R:
        return (r - R)
    else: 
        return 1

def matSci(lmbda,mat):
    #Returns the sigma, albedeo, and cos(theta) grain values based on the wavelength and materials (*no material dependence yet*)
    #Uses linear interpolation to return parameters
    if mat == 'test':
        albedo = 0.6
        g = 0.6
        sigma = 10*0.6 / (1000 * 0.5 * 3.086e18)
        
    else:
        #Loading model data
        data = np.loadtxt('dustdata.txt')
        lam = data[:,0]
        albedogrid = data[:,1]
        ggrid = data[:,2]
        ext = data[:,3]
        absorp = data[:,4]
       
        #Converts absorption coeffiecients to cm^2 per H nucleon
        absorp *= 1.87e-26
        #Calculates scattering coeffiecients
        sigmagrid = ext - absorp
       
        #Interpolation of sigma, albedo, g using model grids
        falb = interpolate.interp1d(lam,albedogrid)
        albedo = float(falb(lmbda))
        
        fg = interpolate.interp1d(lam,ggrid)
        g = float(fg(lmbda))
        
        fsig = interpolate.interp1d(lam,sigmagrid)
        sigma = float(fsig(lmbda))
        
    return sigma,albedo,g

def genKobs(Aobs,L):
    psphere = np.loadtxt('pointsonsphere.txt')
    
    Kobs = np.zeros([psphere.shape[0],2])
    
    for j in range(psphere.shape[0]):
        r = np.sqrt((Aobs[0]-psphere[j,0])**2+(Aobs[1]-psphere[j,1])**2+(Aobs[2]-psphere[j,2])**2)
        Kobs[j,0] = np.arccos((psphere[j,2]-Aobs[2])/(r))
            
        if np.abs(Kobs[j,0]) < 1e-6 or np.abs(Kobs[j,0] - np.pi) < 1e-6:
            Kobs[j,1] = 0.0
#            print(Kobs[j,:])
        else:
            Kobs[j,1] = np.arcsin((psphere[j,1]-Aobs[1])/((r)*np.sin(Kobs[j,0])))
            if psphere[j,0] < L/2:
                Kobs[j,1] = np.pi - Kobs[j,1]
#            print(Kobs[j,:])

    return Kobs

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
    npos[0] = pos[0] + step*np.sin(theta)*np.cos(phi)
    npos[1] = pos[1] + step*np.sin(theta)*np.sin(phi)
    npos[2] = pos[2] + step*np.cos(theta)
    return npos

def interpolateDen(pos,density,axmin,axmax,C,khat):
    if sky(pos,axmax,1e-2) != 0 and sky(pos,axmax,1e-2) !=1:
        return 0

    x,y,z = pos[0],pos[1],pos[2]
    
    i = int((x-axmin)/(axmax-axmin)*(C-1.0))
    j = int((y-axmin)/(axmax-axmin)*(C-1.0))
    k = int((z-axmin)/(axmax-axmin)*(C-1.0))
    
    if sky(pos,axmax,1e-2) == 0:
        return density[i,j,k]
        
    posPrime = posUpdate(pos,khat,axmax/C)
    xprime,yprime,zprime = posPrime[0],posPrime[1],posPrime[2]

    iprime = int((xprime-axmin)/(axmax-axmin)*(C-1.0))
    jprime = int((yprime-axmin)/(axmax-axmin)*(C-1.0))
    kprime = int((zprime-axmin)/(axmax-axmin)*(C-1.0))
    
    if sky(posPrime,axmax,1e-2) != 1:
        return density[i,j,k]
        
    avgDen = (density[i,j,k] + density[iprime,j,k] + density[i,jprime,k] + density[i,j,kprime] + density[iprime,jprime,k] + density[iprime,j,kprime] + density[i,jprime,kprime] + density[iprime,jprime,kprime])/8
    #print(avgDen)
    return np.abs(avgDen)

def odsSample(pos,sigma,omega,g,density,tol,khat,C,L):
    oda = 0.0 
    
    while sky(pos,L,tol):
        p = np.random.rand(1)    
        ods = -np.log(p)
        odsp = 0.0
        
        while np.abs(ods-odsp) > tol:
            den = interpolateDen(pos,density,0.0,L,C,khat)
            if den == 0:
                step = L/C 
            else:
                step = ods / (sigma*den) 
#            print(den)
#            print(step)
            
            poscheck = posUpdate(pos,khat,step)
            if sky(poscheck,L,tol) != 1:
                step -= sky(poscheck,L,tol)
                odsp += sigma * den * step
                oda += (1/omega - 1)*odsp
                pos = posUpdate(pos,khat,step)
#                print('sky',pos)
                return oda
                
            elif step < L/C:
                odsp += sigma * den * step
                pos = posUpdate(pos,khat,step)
#                print('<',pos)
                break
            
            else:    
                odsp += sigma * den * L/C
                pos = posUpdate(pos,khat,L/C)
                ods -= sigma * den * L/C
#                print('>',pos)
            
        oda += (1/omega - 1)*odsp
        khat = genDirection(khat,g)
        
    return oda

def monteCarlo(density,mat,lmbda,Aobs,M,tol,L,C):
    t0 = time.time()
    sigma,omega,g = matSci(lmbda,mat)
    intensities = np.zeros(Aobs.shape[0])        
    
    for a in range(Aobs.shape[0]):
        Kobs = genKobs(Aobs[a,:],L)
        K = Kobs.shape[0]
        W = np.zeros([K,M])
        for k in range(K): 
            for m in range(M):
                pos = Aobs[a,:]
                oda = odsSample(pos,sigma,omega,g,density,tol,Kobs[k],C,L)
                W[k,m] = np.exp(-oda)    
                    
        intensities[a] = 1/(K*M)*np.sum(W)           
    print('time',(time.time()-t0))
    return intensities
