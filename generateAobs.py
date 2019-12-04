## GENERATES RANDOM OBSERVATION POINTS BASED ON DENSITY OF CLOUD ##

#HISTOGRAM RESOLUTION
res = 63
#EXPONENT OF WEIGHT (EMPHASIZES DIFFERENCE BETWEEN HIGH AND LOW WEIGHTS)
weightpower = 4.

import numpy as np
import matplotlib.pyplot as plt

#Converting cloud densities to parsec numbers
cloud = np.load('sampledensity.npy')*1e53/(3.086e18/64)**3

grid = np.arange(64)
x1d = np.zeros(64)
y1d = np.zeros(64)
z1d = np.zeros(64)

# Sum of points along each axis
for i in range(64):
    x1d[i] = np.sum(cloud[i,:,:])
    y1d[i] = np.sum(cloud[:,i,:])
    z1d[i] = np.sum(cloud[:,:,i])
    
# Calculating actual weights
xweights = (x1d/np.sum(x1d))
yweights = (y1d/np.sum(y1d))
zweights = (z1d/np.sum(z1d))
    
#Plotting cloud weights along each axis
fig1 = plt.figure(figsize = (15,4))

plt.subplot(1,3,1)
plt.plot(grid,xweights)
plt.xlabel('x')
plt.ylabel('frequency')

plt.subplot(1,3,2)
plt.plot(grid,yweights)
plt.xlabel('y')

plt.subplot(1,3,3)
plt.plot(grid,zweights)
plt.xlabel('z')
plt.show()

# Calculating powered weights
x1d = x1d**weightpower
y1d = y1d**weightpower
z1d = z1d**weightpower

xweights = (x1d/np.sum(x1d))
yweights = (y1d/np.sum(y1d))
zweights = (z1d/np.sum(z1d))

# Generating random weighted points
randx = np.zeros(10000)
randy = np.zeros(10000)
randz = np.zeros(10000)
iters = 0
while iters < 10000:
    potx = np.random.choice(grid,1,p = xweights)
    poty = np.random.choice(grid,1,p = yweights)
    potz = np.random.choice(grid,1,p = zweights)
    r = np.sqrt((potx-32.)**2. + (poty-32.)**2. + (potz-32)**2.)
    if r <= 32.:
        randx[iters] = potx
        randy[iters] = poty
        randz[iters] = potz
        iters += 1
        
# Creates a 1d array containing (x,y,z) coordinates of each generated point
points = np.zeros((10000,3))
for i in range(10000):
    points[i] = np.array([randx[i], randy[i], randz[i]])
    
#Plotting histograms of generated points along each axis
fig2 = plt.figure(figsize = (15,4))

plt.subplot(1,3,1)
plt.hist(randx,res)
plt.xlabel('x')
plt.ylabel('points generated')

plt.subplot(1,3,2)
plt.hist(randy,res)
plt.xlabel('y')

plt.subplot(1,3,3)
plt.hist(randz,res)
plt.xlabel('z')
plt.show()

#Plotting generated points along each plane
fig3 = plt.figure(figsize = (15,4))

plt.subplot(1,3,1)
plt.plot(randx,randy,'.')
plt.xlabel('x')
plt.ylabel('y')

plt.subplot(1,3,2)
plt.plot(randz,randx,'.')
plt.xlabel('z')
plt.ylabel('x')

plt.subplot(1,3,3)
plt.plot(randy,randz,'.')
plt.xlabel('y')
plt.ylabel('z')

plt.show()

np.savetxt('Aobs.txt',points)
