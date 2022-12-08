import numpy as np
import numba as nb
from particle_class import *
import time
import imageio
from matplotlib import pyplot as plt
from scipy import fft
plt.ion()


parts=particles(npart=1, n = 100, soft=1, periodic=False)

parts.single_particle()

parts.get_kernel()
xy=parts.x.copy()
parts.get_pot()

rho=parts.rho.copy()
pot=parts.pot.copy()
osamp=40


fig = plt.figure()
ax = fig.add_subplot(111)
crap=ax.imshow(parts.rho[:parts.ngrid,:parts.ngrid]**0.5)
times = 500
x_data = np.empty(500)
y_data = np.empty(500)

for i in range(500):

    for j in range(osamp):
        parts.take_step(dt=0.01)

    print("Step", i)
    x_data[i] = parts.x[0, 0]
    y_data[i] = parts.x[0, 1]
    

    crap.set_data(parts.rho[:parts.ngrid,:parts.ngrid]**0.5)
    plt.savefig(f'figs/single_particle/{i:003}', dpi = 50)
    plt.pause(0.001)

plt.clf()
plt.plot(x_data, label = "X Position")
plt.plot(y_data, label = "Y Position")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Position")
plt.savefig("figs/single_position.jpg")
plt.show()
