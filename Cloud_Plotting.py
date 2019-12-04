import argparse                  # allows us to deal with arguments to main()
from argparse import RawTextHelpFormatter
from Cloud import Cloud
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
from MonteCarloScheme import monteCarlo
import multiprocessing as mp
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("D",type=float,
                        help="Fractal Dimention")
    parser.add_argument("L",type=float,
                        help="Size of cube axis")
    parser.add_argument("C",type=int,
                        help="Number of Cubes Axis")
    parser.add_argument("N",type=int,
                        help="Generation particle number")
    parser.add_argument("mat",type=str,
                        help="Material of Cloud")
    parser.add_argument("lmbd",type=float,
                        help="Wavelength")
    parser.add_argument("uniform",type=int,
                        help="Uses uniform density if True (1)")
    parser.add_argument("sample",type=int,
                        help="Uses sample density if True (1)")
    
    args = parser.parse_args()
    D = args.D
    L = args.L
    C = args.C
    N = args.N
    mat = args.mat
    lmbd = args.lmbd 
    uniform = args.uniform
    sample = args.sample

    cloud, density = Cloud(D,N,L,C,uniform = uniform, sample = sample).point_array
    
#    print(density)
    #Plotting For the Particles in Space
     
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#
#    ax.scatter(cloud[:,0], cloud[:,1], cloud[:,2], c='r', marker='.')
#    r = [0, L]
#    for s, e in combinations(np.array(list(product(r, r, r))), 2):
#        if np.sum(np.abs(s-e)) == r[1]-r[0]:
#            ax.plot3D(*zip(s, e), color="b")
#        
#    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
#    x = L/2*np.cos(u)*np.sin(v) + L/2
#    y = L/2*np.sin(u)*np.sin(v) + L/2
#    z = L/2*np.cos(v) + L/2
#    ax.plot_wireframe(x, y, z, color="gray",alpha = 0.5)
#
#    plt.show()
    
    R = L/2
    if mat == 'test':
        Aobs = np.zeros([20,3]) + R
        for n in range(20):
            Aobs[n,2] += 0.6*0.5*n / ((10*0.6 / (1000 * 0.5 * 3.086e18))*1000)
    else:
        Aobs = np.loadtxt('Aobs.txt')
        Aobs *= (L/C)
        
    Aobs_size = Aobs.shape[0]
    
#    print(Aobs)
    Obs_list = []
    pross = mp.cpu_count()    
    Obs_list = np.array_split(Aobs,pross)

#    print(Obs_list)
    
    tasks = []
    for i in range(pross):
        tuple = (density,mat,lmbd,Obs_list[i],10,1e-2,L)
        tasks.append(tuple)
    
#    print(tasks)
    
    myPool = mp.Pool(pross)
    intensities = myPool.starmap(monteCarlo, tasks)
    
    intensityArr = np.concatenate(intensities).ravel()
    np.savetxt('testintensity.txt',intensityArr)
    
#    intensity = monteCarlo(density,mat,lmbd,Aobs,10,1e-2,L)

#    np.savetxt('testintensity.txt',intensity)
