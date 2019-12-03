import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

def h_l(l,g,omega):
    sig = g
    if l is 0:
        sig = 1

    return((2*float(l)+1)*(1-omega*sig))

def p_l(k,h,l):
    

g = 0.4
omega = 0.5
# omega = np.linspace(1e-6,1,endpoint=True,num=1e4)

L = [0,1,2,3,4,5]

sumh = 0
for l in L:
    sumh += funch(l,g,omega)

print(sumh)



