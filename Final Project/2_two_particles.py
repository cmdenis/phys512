import numpy as np
import numba as nb
from particle_class import *
from matplotlib import pyplot as plt
from scipy import fft
plt.ion()



# Initializing our system with 2 particles
parts=Particles(npart=2, n = 100, soft=2 ,periodic=True)
parts.two_particles()

# Get the kernel
parts.get_kernel()

# Copy 
xy=parts.x.copy()
parts.get_pot()
rho=parts.rho.copy()
pot=parts.pot.copy()

osamp=400    


fig = plt.figure()
ax = fig.add_subplot(111)
res=ax.imshow(parts.rho[:parts.ngrid,:parts.ngrid]**0.5)

times = 500

energies = np.empty([3, times])
print(parts.m)

for i in range(times):

    for j in range(osamp):
        parts.take_step(dt=0.001)
        res.set_data(parts.rho[:parts.ngrid, :parts.ngrid])


    # Energies
    kin = np.sum(parts.v**2 * parts.m)/2
    pot = np.sum(parts.rho * parts.pot)
    tot = kin+pot
    # Storing values in array for plotting
    energies[0, i] = kin 
    energies[1, i] = pot 
    energies[2, i] = tot 
    print("\nKinetic Energy:", kin, "\nPotential:", pot, "\nTotal Energy:", tot)


    res.set_data(parts.rho[:parts.ngrid, :parts.ngrid])
    plt.savefig(f'figs/two_particles/{i:003}', dpi = 50)
    plt.pause(0.001)

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
