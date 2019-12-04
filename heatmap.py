###PLOTTING OPTIONS
color = 'Greys'
resolution = 64
###

import numpy as np
import matplotlib.pyplot as plt

#Converting cloud densities to parsec numbers
cloud = np.load('sampledensity.npy')*1e53/(3.086e18/64)**3

#Setting up 64x64 grid for cloud
x = np.arange(64)
y = np.arange(64)

X,Y = np.meshgrid(x,y)

#Summing the cloud densities of each axis
xsum = np.zeros((64,64))
ysum = np.zeros((64,64))
zsum = np.zeros((64,64))
for i in range(64):
    for j in range(64):
        xsum[i,j] = np.sum(cloud[:,i,j])
        ysum[i,j] = np.sum(cloud[j,:,i])
        zsum[i,j] = np.sum(cloud[i,j,:])

#Circular cloud boundary for plotting
theta = np.linspace(0,2*np.pi,1000)
xcirc = 31*np.cos(theta) + 31.5
ycirc = 31*np.sin(theta) + 31.5

#Plotting each density plane
fig = plt.figure(figsize = (15,4))

#XY plane
plt.subplot(1,3,1)
plt.contourf(X,Y,zsum,resolution,cmap=color)
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.plot(xcirc,ycirc,'r')

#XZ plane
plt.subplot(1,3,2)
plt.contourf(X,Y,ysum,resolution,cmap=color)
plt.colorbar()
plt.xlabel('z')
plt.ylabel('x')
plt.plot(xcirc,ycirc,'r')

#YZ plane
plt.subplot(1,3,3)
plt.contourf(X,Y,xsum,resolution,cmap=color)
plt.colorbar()
plt.xlabel('y')
plt.ylabel('z')
plt.plot(xcirc,ycirc,'r')

plt.tight_layout()
