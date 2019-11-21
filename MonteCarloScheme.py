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
        return (r-R)
    else: 
        return 1

def matSci(lmbda,mat):
    #Returns the sigma, albedeo, and cos(theta) grain values based on the wavelength and materials (*no material dependence yet*)
    #Uses linear interpolation to return parameters
    if mat == 'test':
        albedo = 0.6
        g = 0.6
        sigma = 4.1666
        
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
        #Calculates scattering coeeffiecients
        sigmagrid = ext - absorp
       
        #Interpolation of sigma, albedo, g using model grids
        falb = interpolate.interp1d(lam,albedogrid)
        albedo = falb(lmbda)
        
        fg = interpolate.interp1d(lam,ggrid)
        g = fg(lmbda)
        
        fsig = interpolate.interp1d(lam,sigmagrid)
        sigma = fsig(lmbda)
        
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
    npos[0] = pos[0] + step*np.cos(theta)*np.sin(phi)
    npos[1] = pos[1] + step*np.sin(theta)*np.sin(phi)
    npos[2] = pos[2] + step*np.cos(theta)
    return npos

def interpolateDen(pos,density,axmin,axmax,C,khat):
    x,y,z = pos[0],pos[1],pos[2]
    
    i = int((x-axmin)/(axmax-axmin)*(C-1.0))
    j = int((y-axmin)/(axmax-axmin)*(C-1.0))
    k = int((z-axmin)/(axmax-axmin)*(C-1.0))
    
    posPrime = posUpdate(pos,khat,np.sqrt(3))
    xprime,yprime,zprime = posPrime[0],posPrime[1],posPrime[2]
    
    iprime = int((xprime-axmin)/(axmax-axmin)*(C-1.0))
    jprime = int((yprime-axmin)/(axmax-axmin)*(C-1.0))
    kprime = int((zprime-axmin)/(axmax-axmin)*(C-1.0))
    
    if iprime >= 64:
        iprime = i 
    if jprime >= 64:
        jprime = j 
    if kprime >= 64:
        kprime = k
    
#    avgDen = density[i,j,k]*(1-Wx)*(1-Wy)*(1-Wz) + density[iprime,j,k]*(Wx)*(1-Wy)*(1-Wz) 
#    + density[i,jprime,k]*(1-Wx)*(Wy)*(1-Wz) + density[i,j,kprime]*(1-Wx)*(1-Wy)*(Wz) 
#    + density[iprime,jprime,k]*(Wx)*(Wy)*(1-Wz) + density[iprime,j,kprime]*(Wx)*(1-Wy)*(Wz)
#    + density[i,jprime,kprime]*(1-Wx)*(Wy)*(Wz) + density[iprime,jprime,kprime]*(Wx)*(Wy)*(Wz)
        
    avgDen = density[i,j,k]*0.79 + density[iprime,j,k]*0.03 + density[i,jprime,k]*0.03 + density[i,j,kprime]*0.03 + density[iprime,jprime,k]*0.03 + density[iprime,j,kprime]*0.03 + density[i,jprime,kprime]*0.03 + density[iprime,jprime,kprime]*0.03
    #print(avgDen)
    return np.abs(avgDen)

def odsSample(pos,sigma,omega,g,density,tol,khat):
    oda = 0.0 
    
    C = 64.0
    L = 64.0
    it = 0 
    
    while sky(pos,L,tol):
        p = np.random.rand(1)    
        ods = -np.log(p)
        odsp = 0.0
        
        while np.abs(ods-odsp) > tol:
            den = interpolateDen(pos,density,0.0,L,C,khat)
            if den == 0:
                step = np.sqrt(3)
            else:
                step = ods / (sigma*den) 
#            print(den)
#            print(step)
            
            poscheck = posUpdate(pos,khat,step)
            if sky(poscheck,L,tol) != 0 and sky(poscheck,L,tol) !=1:
                step -= sky(poscheck,L,tol)
                odsp += sigma * den * step
                pos = posUpdate(pos,khat,step)
#                print(pos)
                break
                
            elif step < np.sqrt(3):
                odsp += sigma * den * step
                pos = posUpdate(pos,khat,step)
#                print(pos)
                break
            
            else:    
                odsp += sigma * den * np.sqrt(3)
                pos = posUpdate(pos,khat,np.sqrt(3))
                ods -= sigma * den * np.sqrt(3)
#                print(pos)
            
        oda += (1/omega - 1)*odsp
        khat = genDirection(khat,g)
        it+=1
        #print(it)
        
    return oda

def monteCarlo(density,mat,lmbda,Aobs,Kobs,M,tol):
    t0 = time.time()
    sigma,omega,g = matSci(lmbda,mat)
    K = Kobs.shape[0]
    intensities = np.zeros(Aobs.shape[0])        
    
    for a in range(Aobs.shape[0]):
        W = np.zeros([K,M])
        for k in range(K): 
            for m in range(M):
                pos = Aobs[a,:]
                oda = odsSample(pos,sigma,omega,g,density,tol,Kobs[k])
                W[k,m] = oda    
                    
        intensities[a] = 1/(K*M)*np.sum(W)           
    print('time',(time.time()-t0))
    return intensities

density = np.zeros([64,64,64]) + 0.15
Aobs = np.zeros([1,3]) + 64.0/2
#print(Aobs)
M = 1
kobs = np.array([[0,0]])#,
#                 [np.pi/2,0],
#                 [np.pi/2,np.pi/2],
#                 [np.pi/2,np.pi],
#                 [np.pi/2,-np.pi/2],
#                 [np.pi/2,-np.pi]])

print('intensitiy' , monteCarlo(density,'test',400e-5,Aobs,kobs,M,1e-2))
