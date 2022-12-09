import numpy as np
import numba as nb
from particle_class import *
from matplotlib import pyplot as plt
from scipy import fft
plt.ion()




# We simulate the periodic boundary conditions first

# Instantiate our system
npart = 200000  # Number of particles
ng = 1000       # Size of grid
parts = Particles(
    npart=npart, 
    n = ng, 
    periodic = True
    )

# Initialize particles
parts.many_particles()
# Get the kernel
parts.get_kernel()
# Copy stuff to avoid problems
xy = parts.x.copy()
parts.get_pot()
rho = parts.rho.copy()
pot = parts.pot.copy()

# Over sampling
over_sample = 12    

# Setup figure
fig = plt.figure()
ax = fig.add_subplot(111)
res = ax.imshow(parts.rho[:parts.ngrid, :parts.ngrid]**0.5)

# Number of frames
times = 500

# Energy arrays
energies_per = np.empty([3, times])
energies_nonper = np.empty([3, times])





# Main loop
for i in range(times):

    # Oversampling loop
    for j in range(over_sample):
        parts.take_step(dt = 0.01)

    # Energies
    kin = np.sum(parts.v**2)/2
    pot = np.sum(parts.rho * parts.pot)
    tot = kin+pot

    # Storing values in array for plotting
    energies_per[0, i] = kin 
    energies_per[1, i] = pot 
    energies_per[2, i] = tot 
    print("\nKinetic Energy:", kin, "\nPotential:", pot, "\nTotal Energy:", tot)


    res.set_data(parts.rho[:parts.ngrid, :parts.ngrid])
    plt.savefig(f'figs/periodic/{i:003}', dpi = 50)
    plt.pause(0.001)

















# We simulate the same system but with nonperiodic conditions

# Initializing our system with 400 particles
parts = Particles(
    npart=npart, 
    n = ng, 
    periodic = False
    )

# Initialize particles
parts.many_particles()
# Get the kernel
parts.get_kernel()
# Get potential
parts.get_pot()
# Copy stuff to avoid problems
xy = parts.x.copy()
rho = parts.rho.copy()
pot = parts.pot.copy()

# Setting plotting stuff
plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111)
res = ax.imshow(parts.rho[:parts.ngrid,:parts.ngrid]**0.5)




# Main loop
for i in range(times):

    # Oversampling loop
    for j in range(over_sample):
        parts.take_step(dt=0.01)

    # Energies
    kin = np.sum(parts.v**2)/2
    pot = np.sum(parts.rho * parts.pot)
    tot = kin+pot
    # Storing values in array for plotting
    energies_nonper[0, i] = kin 
    energies_nonper[1, i] = pot 
    energies_nonper[2, i] = tot 
    print("\nKinetic Energy:", kin, "\nPotential:", pot, "\nTotal Energy:", tot)

    # Plotting and saving
    res.set_data(parts.rho[:parts.ngrid, :parts.ngrid])
    plt.savefig(f'figs/non_periodic/{i:003}', dpi = 50)
    plt.pause(0.001)









# Plotting the evolution of total energy for periodic AND non-periodic systems
plt.clf()
plt.plot(energies_per[2, :], label = "Total Energy (Periodic)")
plt.plot(energies_nonper[2, :], label = "Total Energy (Non-Periodic)")
plt.xlabel("Timesteps")
plt.ylabel("Energy")
plt.legend()
plt.title("Comparison of Energy For Periodic and Non-Periodic BC")
plt.savefig("figs/per_vs_nonper.jpg")
plt.show()