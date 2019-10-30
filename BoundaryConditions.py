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
    
    boundx = np.linspace(0.0,L,C+1)
    boundy = np.linspace(0.0,L,C+1)
    boundz = np.linspace(0.0,L,C+1)
    print(boundx)
    
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

def main():
    
    arr = np.array([[0.5,0.5,0.5],
                    [1.5,0.5,0.5],
                    [0.5,1.5,0.5],
                    [0.5,0.5,1.5],
                    [1.5,1.5,0.5],
                    [1.5,0.5,1.5],
                    [0.5,1.5,1.5],
                    [1.5,1.5,1.5]])
    
    #print(arr[:,0],arr[:,1],arr[:,2])
    print(densityFunc(np.squeeze(arr[:,0]),np.squeeze(arr[:,1]),np.squeeze(arr[:,2]),2,2))
    
main()
