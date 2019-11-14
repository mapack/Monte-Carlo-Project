from Cloud import Cloud
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations

D = 2.6
N = 5
L = 64
C = 64

cloudR = Cloud(D,N,L,C).point_array
cloud = cloudR[0]
print('Density Array: ' + str(cloudR[1][0]))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(cloud[:,0], cloud[:,1], cloud[:,2], c='r', marker='.')
r = [0, L]
for s, e in combinations(np.array(list(product(r, r, r))), 2):
    if np.sum(np.abs(s-e)) == r[1]-r[0]:
        ax.plot3D(*zip(s, e), color="b")
        
#u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
#x = L/2*np.cos(u)*np.sin(v) + L/2
#y = L/2*np.sin(u)*np.sin(v) + L/2
#z = L/2*np.cos(v) + L/2
#ax.plot_wireframe(x, y, z, color="gray",alpha = 0.5)

plt.show()
