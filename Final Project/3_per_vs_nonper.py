import numpy as np
import numba as nb
import time
import imageio
from matplotlib import pyplot as plt
from scipy import fft
plt.ion()


def inbound_array_np(xy,n):
    '''
    Applies periodic conditions with numpy arrays.
    If the position goes out of bound it gets put back at the begining.
    '''
    xy[xy<-0.5] = xy[xy <- 0.5] + n
    xy[xy>=n-0.5] = xy[xy >= n-0.5]-n

@nb.njit
def mask_array_np(xy,m,n):
    '''
    Mask for non-periodic boundary conditions.
    If a value goes out of bound, its mass gets reduced to 0.
    '''
    for i in range(xy.shape[1]):
        m[xy[:,i]<-0.5]=0
        m[xy[:,i]>=n-0.5]=0

@nb.njit(parallel=True)
def get_grad(xy,pot,grad, m):
    '''
    Function to compute the gradient of an array.
    It basically uses the slope between neighboring cell sites to create the gradient.
    Beforehand it takes care of some conditionals involving the position of the particles
    with respect to the first and last cell in the density (rho) array.
    '''
    n=pot.shape[0]
    for i in nb.prange(xy.shape[0]):
        if m[i] > 0.0:
            if xy[i,0]<0:
                ix0=n-1
                ix1=0
                fx=xy[i,0]+1
            else:
                ix0=int(xy[i,0])
                ix1=ix0+1
                fx=xy[i,0]-ix0
                if ix1==n:
                    ix1=0
            if xy[i,1]<0:
                iy0=n-1
                iy1=0
                fy=xy[i,1]+1
            else:
                iy0=int(xy[i,1])
                iy1=iy0+1
                fy=xy[i,1]-iy0
                if iy1==n:
                    iy1=0

            grad[i,0] = (pot[ix1,iy1]-pot[ix0,iy1])*fy+(pot[ix1,iy0]-pot[ix0,iy0])*(1-fy)
            grad[i,1] = (pot[ix1,iy1]-pot[ix1,iy0])*fx+(pot[ix0,iy1]-pot[ix0,iy0])*(1-fx)
        else:
            grad[i,0] = 0
            grad[i,1] = 0

@nb.njit
def hist2d_wmass(xy,mat,m):
    '''
    Creates the histogram for the density (rho) array, 
    based on the masses and positions of the particles.
    '''
    nx=xy.shape[0]
    for i in range(nx):
        ix = int(xy[i,0]+0.5)
        iy = int(xy[i,1]+0.5)
        if m[i] > 0: #we can set m=0 to flag particles
            mat[ix, iy] = mat[ix, iy] + m[i]

class particles:
    def __init__(self,npart=10000,n=1000,soft=1,periodic=True):
        self.x = np.empty([npart,2])        # Positions
        self.v = np.empty([npart,2])        # Velocities
        self.f = np.empty([npart,2])        # Forces
        self.m = np.empty(npart)            # Masses
        self.kernel = []                    # Kernel
        self.kernelft = []                  # FFT of kernel
        self.npart = npart                  # Number of particles
        self.ngrid = n                      # Width and height of grid

        if periodic: # Creating array for rho and potential based on the size of the grid
            self.rho = np.empty([self.ngrid,self.ngrid])
            self.pot = np.empty([self.ngrid,self.ngrid])

        else:   # Making rho twice the size if not periodic
            self.rho = np.empty([2*self.ngrid,2*self.ngrid])
            self.pot = np.empty([2*self.ngrid,2*self.ngrid])

        self.soft = soft            # Softening
        self.periodic = periodic    # Periodicity boolean

    def many_particles(self):
        '''
        Initializes the positions of n particles uniformly distributed.
        '''
        np.random.seed(514) # Use a seed for repeatability
        self.x = np.random.uniform(0, self.ngrid, [self.npart, 2])
        self.v[:] = 0    # 0 Initial velocity (Particle 1)
        self.m[:] = 1

    def get_kernel(self):
        '''Create the kernel (and fft of kernel) to convolve with density grid.'''
        ngrid = self.ngrid
        soft = self.soft
        
        if self.periodic:
            x=np.fft.fftfreq(ngrid)*ngrid
            rsqr=np.outer(np.ones(ngrid),x**2)
            rsqr=rsqr+rsqr.T
            rsqr[rsqr<soft**2]=soft**2
            self.kernel=rsqr**-0.5
        else:
            # Makes the kernel twice as big for non-periodic BC
            x=np.fft.fftfreq(ngrid*2)*ngrid*2
            rsqr=np.outer(np.ones(ngrid*2),x**2)
            rsqr=rsqr+rsqr.T
            rsqr[rsqr<soft**2]=soft**2
            self.kernel=rsqr**-0.5

        self.kernelft=fft.rfft2(self.kernel)    # Fourier transform of kernel

    def get_rho(self):
        '''Method to get the density in a grid of the particles.'''
        if self.periodic:
            # If periodic, give appropriate location to particles that are outside bounds
            inbound_array_np(self.x,self.ngrid)
        else:
            # If not periodic, sets masses to 0, if they're outside the bounds
            mask_array_np(self.x, self.m, self.ngrid)
        self.rho[:]=0

        # Creates histogram for density array (rho)
        hist2d_wmass(self.x,self.rho,self.m)
    
    def get_pot(self):
        '''Method to get the potential using a convolution of kernel and density'''
        self.get_rho()
        rhoft=fft.rfft2(self.rho)
        n=self.ngrid
        if not(self.periodic):
            n=n*2
        self.pot=fft.irfft2(rhoft*self.kernelft,[n,n])
    
    def get_forces(self):
        '''Get the force on every particle using the gradient of the potential.'''
        get_grad(self.x, self.pot, self.f, self.m)
    
    def take_step(self,dt=1):
        # Takes step using leapfrog
        self.x[:]=self.x[:]+dt*self.v
        self.get_pot()
        self.get_forces()
        self.v[:]=self.v[:]+self.f*dt


# We simulate the periodic boundary conditions first

# Initializing our system with 400 particles
npart = 700
ng = 200
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

osamp=60    

fig = plt.figure()
ax = fig.add_subplot(111)
res=ax.imshow(parts.rho[:parts.ngrid,:parts.ngrid]**0.5)

times = 500

energies_per = np.empty([3, times])
energies_nonper = np.empty([3, times])
print("Ready to start")
for i in range(times):

    for j in range(osamp):
        parts.take_step(dt=0.01)

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


    res.set_data(parts.rho[:parts.ngrid, :parts.ngrid])
    plt.savefig(f'figs/non_periodic/{i:003}', dpi = 50)
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
plt.savefig("figs/per_vs_nonper.jpg")
plt.show()