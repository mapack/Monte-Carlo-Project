import argparse                  # allows us to deal with arguments to main()
from argparse import RawTextHelpFormatter
from Cloud import Cloud
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
from MonteCarloScheme import monteCarlo

def main():
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
    parser.add_argument("uniform",type=bool,
                        help="Uses uniform density if True")
    
    args = parser.parse_args()
    D = args.D
    L = args.L
    C = args.C
    N = args.N
    mat = args.mat
    lmbd = args.lmbd 
    uniform = args.uniform

    cloud, density = Cloud(D,N,L,C,uniform = uniform).point_array

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
            Aobs[i,2] += i/(49.342*0.019)
    else:
        Aobs = np.zeros([10,000,3])
        
    
    intensity = monteCarlo(density,mat,lmbd,Aobs,10,1e-2)
