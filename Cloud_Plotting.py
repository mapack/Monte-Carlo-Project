from Cloud import Cloud
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

D = 2.6
N = 5
L = 100

cloud = Cloud(D,N,L).point_array

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(cloud[:,0], cloud[:,1], cloud[:,2], c='r', marker='o')
plt.show()
