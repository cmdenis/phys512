import numpy as np
import matplotlib.pyplot as plt
import ctypes
import time
import random

# Test it out with text file
rand_val = np.loadtxt('rand_points.txt')

a = 1.236
b = -0.607

x_data = a*rand_val[:, 0] + b*rand_val[:, 1]
y_data = rand_val[:, 2]

plt.scatter(x_data, y_data, s = 1)
plt.title("Cross-Section of random number plotted in 3D")
plt.savefig("figs/a7q1_cross_section.jpg")
#plt.show()
plt.clf()







# Looking at the same thing but for python's PRNG
n = 300
rand_val_x = np.random.rand(n)
rand_val_y = np.random.rand(n)
rand_val_z = np.random.rand(n)


fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')

ax.scatter(rand_val_x, rand_val_y, rand_val_z)
plt.savefig("figs/a7q1_numpy.jpg")
#plt.show()
plt.clf()










# Get random numbers from local machine
n=300000
mylib=ctypes.cdll.LoadLibrary("libc.dylib")
rand=mylib.rand
rand.argtypes=[]
rand.restype=ctypes.c_int

def get_rands_nb(vals):
    n=len(vals)
    for i in range(n):
        vals[i]=rand()
    return vals

def get_rands(n):
    vec=np.empty(n,dtype='int32')
    get_rands_nb(vec)
    return vec


vec=get_rands(n*3)

vv=np.reshape(vec,[n,3])
vmax=np.max(vv,axis=1)

maxval=5e8

vv2=vv[vmax<maxval, :]

np.savetxt("rand_points_local.txt", vv2) # Saving to machine

rand_val = np.loadtxt('rand_points_local.txt') # Loading them back

xx = rand_val[:, 0]
yy = rand_val[:, 1]
zz = rand_val[:, 2]


a = 1.236
b = -0.607

x_data = a*rand_val[:, 0] + b*rand_val[:, 1]
y_data = rand_val[:, 2]

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')

ax.scatter(xx, yy, zz, s = 1)
#plt.savefig("figs/a7q1_numpy.jpg")
plt.show()
plt.clf()

plt.scatter(x_data, y_data, s = 1)
plt.title("Cross-Section of random number plotted in 3D")
plt.savefig("figs/a7q1_cross_section_local.jpg")
plt.show()
plt.clf()



