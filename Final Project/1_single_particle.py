import numpy as np
import numba as nb
from particle_class import *
import time
import imageio
from matplotlib import pyplot as plt
from scipy import fft
plt.ion()

# Instantiating the system
parts = Particles(
    npart = 1,          # Number of particles
    n = 100,            # Size of grid
    soft = 1,           # Softening
    periodic = False    # Periodic Boolean
    )

# Initializing the system
parts.single_particle()

# Getting the kernel
parts.get_kernel()

# Getting the potential
parts.get_pot()

# Copying stuff to avoid problems
xy = parts.x.copy()
rho = parts.rho.copy()
pot = parts.pot.copy()


# Over sampling
over_sample = 12  

# Setting figure stuff
fig = plt.figure()
ax = fig.add_subplot(111)
res = ax.imshow(parts.rho[:parts.ngrid,:parts.ngrid]**0.5)

times = 500                 # Number of frames
x_data = np.empty(times)    # Creating array for x-position
y_data = np.empty(times)    # Creating array for y-position

# Main loop for frames
for i in range(times):

    # Oversampling loop
    for j in range(over_sample):
        parts.take_step(dt = 0.01)

    print("Step", i)
    # Storing positions
    x_data[i] = parts.x[0, 0]
    y_data[i] = parts.x[0, 1]
    
    # Plotting and saving
    res.set_data(parts.rho[:parts.ngrid,:parts.ngrid]**0.5)
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
