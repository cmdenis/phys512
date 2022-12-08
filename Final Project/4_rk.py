import numpy as np
import numba as nb
import time
from particle_class import *
import imageio
from matplotlib import pyplot as plt
from scipy import fft
plt.ion()

# We simulate the periodic boundary conditions first

# Initializing our system
npart = 200000
ng = 1000
parts=particles(npart=npart, n = ng, periodic=True)
parts.many_particles()
print("System Initialized")
# Get the kernel
parts.get_kernel()
print("Got kernel")
# Copy 
xy=parts.x.copy()
parts.get_pot()
rho=parts.rho.copy()
pot=parts.pot.copy()

print(parts.m)

osamp = 3    

fig = plt.figure()
ax = fig.add_subplot(111)
res=ax.imshow(parts.rho[:parts.ngrid,:parts.ngrid]**0.5)

times = 500

energies_per = np.empty([3, times])
energies_nonper = np.empty([3, times])
print("Ready to start")
for i in range(times):

    for j in range(osamp):
        parts.take_step_rk(dt=0.04)

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
    plt.savefig(f'figs/rk_periodic/{i:003}', dpi = 50)
    plt.pause(0.001)


# We simulate the same system but with nonperiodic conditions

# Initializing our system with 400 particles
parts=particles(npart=npart, n = ng, periodic=False)
parts.many_particles()

# Get the kernel
parts.get_kernel()

# Copy 
xy=parts.x.copy()
parts.get_pot()
rho=parts.rho.copy()
pot=parts.pot.copy()


plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111)
res=ax.imshow(parts.rho[:parts.ngrid,:parts.ngrid]**0.5)

for i in range(times):

    for j in range(osamp):
        parts.take_step_rk(dt=0.04)

    # Energies
    kin = np.sum(parts.v**2)/2
    pot = np.sum(parts.rho * parts.pot)
    tot = kin+pot
    # Storing values in array for plotting
    energies_nonper[0, i] = kin 
    energies_nonper[1, i] = pot 
    energies_nonper[2, i] = tot 
    print("\nKinetic Energy:", kin, "\nPotential:", pot, "\nTotal Energy:", tot)


    res.set_data(parts.rho[:parts.ngrid, :parts.ngrid])
    plt.savefig(f'figs/rk_nonperiodic/{i:003}', dpi = 50)
    plt.pause(0.001)


print("here")


plt.clf()
#plt.plot(energies_per[0, :], label = "Kinetic Energy (Periodic)")
#plt.plot(energies[1, :], label = "Potential Energy")
plt.plot(energies_per[2, :], label = "Total Energy (Periodic)")
plt.plot(energies_nonper[2, :], label = "Total Energy (Non-Periodic)")
plt.xlabel("Timesteps")
plt.ylabel("Energy")
plt.legend()
plt.title("Comparison of Energy For Periodic and Non-Periodic BC")
plt.savefig("figs/rk_per_vs_nonper.jpg")
plt.show()