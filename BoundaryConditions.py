import numpy as np

#Tests if each hard wall boundary of the square has been broken.
#Uses the array generated in each generation step after the first
#Requires the length of the cube.
def hardWallBound(arr,L):
    x = arr[:,0]
    y = arr[:,1]
    z = arr[:,3]
    
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
            
        return x,y,z
        
def densityFunc(x,y,z,L,C):
    #C is the number of cubes per axis length L
    density = np.zeros([C,C,C])
    
    for ...:
        zcond
        for...:
            xcond
            for...:
                ycond
                
                if np.all(...):
                    density[...] = np.sum(...)
    
    return density 