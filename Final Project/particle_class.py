import numpy as np
import numba as nb
import time
import imageio
from matplotlib import pyplot as plt
from scipy import fft
#plt.ion()


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

class Particles:
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

    def single_particle(self):
        '''
        Initializes the position of 1 particle at roughly the center of the grid
        '''
        self.x[:] = [self.ngrid//2, self.ngrid//2]
        self.v[:] = 0    
        self.m[:] = 1

    def two_particles(self):
        '''
        Initializes the positions of two particles at the center of the grid 
        with some opposite velocities.
        '''
        n = self.ngrid  # Size of grid

        # Positioning the particles at the center of the plane separated from each other
        self.x[0] = np.array([n/2, n*0.4])  # Positioned at 40% of the grid (Particle 1)
        self.x[1] = np.array([n/2, n*0.6])  # Positioned at 60% of the grid (Particle 2)
        factor = 0.5
        self.v[0] = np.array([-1, 0])*factor    # Initial velocity (Particle 1)
        self.v[1] = np.array([1, 0])*factor     # Initial velocity (Particle 2)

        self.m[:] = 20

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

    def take_step_rk(self, dt=1):
        # Takes step using Runge-Kutta
        # using the second order derivative description here: https://fr.wikipedia.org/wiki/M%C3%A9thodes_de_Runge-Kutta

        initial_x = self.x.copy()

        # finding k1
        self.get_pot()
        self.get_forces()
        k1 = self.f.copy()

        # Finding k2
        self.x[:] = self.x[:] + dt*self.v/2 
        self.get_pot()
        self.get_forces()
        k2 = self.f.copy()

        # Finding k3
        self.x[:] = self.x[:] + dt**2 * k1/4
        self.get_pot()
        self.get_forces()
        k3 = self.f.copy()

        # Finding k4
        self.x[:]=self.x[:] + dt*self.v/2 + dt**2 * (k2/2 - k1/4)
        self.get_pot()
        self.get_forces()
        k4 = self.f.copy()

        self.f = k1

        self.x[:] = initial_x + dt*self.v + (dt**2)*(k1 + k2 + k3)/6

        self.v[:] = self.v[:] + dt*(k1 + 2*k2 + 2*k3 + k4)/6

