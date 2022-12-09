import numpy as np
import numba as nb
from particle_class import *
from matplotlib import pyplot as plt
from scipy import fft
plt.ion()



# Instantiating our system with 2 particles
parts = Particles(
    npart = 2,      # Number of particles
    n = 100,        # Size of grid
    soft = 2,       # Softening
    periodic = True # Periodicity
    )

# Initializing our system
parts.two_particles()

# Get the kernel
parts.get_kernel()

# Get potential
parts.get_pot()

# Copy stuff to avoid problems
xy = parts.x.copy()
rho = parts.rho.copy()
pot = parts.pot.copy()

# Oversampling
over_sample = 400    

# Setting figure up
fig = plt.figure()
ax = fig.add_subplot(111)
res = ax.imshow(parts.rho[:parts.ngrid,:parts.ngrid]**0.5)

times = 500     # Number of frames

# Setting energy array
energies = np.empty([3, times])









# Main loop for frames
for i in range(times):

    # Oversampling loop
    for j in range(over_sample):
        parts.take_step(dt = 0.001)

    # Energies
    kin = np.sum(parts.v**2 * parts.m)/2
    pot = np.sum(parts.rho * parts.pot)
    tot = kin+pot

    # Storing values in array for plotting
    energies[0, i] = kin 
    energies[1, i] = pot 
    energies[2, i] = tot 
    print("\nKinetic Energy:", kin, "\nPotential:", pot, "\nTotal Energy:", tot)

    # Potting and saving frames
    res.set_data(parts.rho[:parts.ngrid, :parts.ngrid])
    plt.savefig(f'figs/two_particles/{i:003}', dpi = 50)
    plt.pause(0.001)








# Plotting stuff for energies over time
plt.clf()
plt.plot(energies[0, :], label = "Kinetic Energy")
plt.plot(energies[1, :], label = "Potential Energy")
plt.plot(energies[2, :], label = "Total Energy")
plt.xlabel("Timesteps")
plt.ylabel("Energy")
plt.legend()
plt.title("Evolution of energy of the system for two particles")
plt.savefig("figs/2particles_energy.jpg")
plt.show()
