import numpy as np
import numba as nb
import time
from matplotlib import pyplot as plt
from scipy import fft



@nb.njit
def rho_maker(xy, m, mat):
    '''Compute the density matrix'''
    nx=xy.shape[0] # Number of particles
    for i in range(nx):
        ix=int(xy[i,0]+0.5)
        iy=int(xy[i,1]+0.5)
        #if m[i]>0: #we can set m=0 to flag particles
        mat[ix, iy] = mat[ix, iy] + m[i] # Filling up the density matrix
    


@nb.jit(parallel = True)
def get_grad1(pos, pot, grad):
    n = pot.shape[0]
    for i in nb.prange(pos.shape[0]):
        # First we do the x axis
        x0 = int(pos[i, 0])
        x1 = np.mod(x0 + 1, n)
        dx0 = pos[i, 0] - x0

        # Second we do the y axis
        y0 = int(pos[i, 1])
        y1 = np.mod(y0 + 1, n)
        dy0 = pos[i, 0] - y0

        # Now we find the gradient
        grad[i, 0] = (pot[x1, y1] - pot[x0, y1]) * dy0 + (pot[x1, y0] - pot[x0, y0]) * (1 - dy0)
        grad[i, 1] = (pot[x1, y1] - pot[x1, y0]) * dx0 + (pot[x0, y1] - pot[x0, y0]) * (1 - dx0)

@nb.njit(parallel=True)
def get_grad(xy,pot,grad):
    n=pot.shape[0]
    for i in nb.prange(xy.shape[0]):
        if xy[i,0]<0:
            print("here:", xy[i,:])
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
        
        grad[i,0]=(pot[ix1,iy1]-pot[ix0,iy1])*fy+(pot[ix1,iy0]-pot[ix0,iy0])*(1-fy)
        grad[i,1]=(pot[ix1,iy1]-pot[ix1,iy0])*fx+(pot[ix0,iy1]-pot[ix0,iy0])*(1-fx)



def make_kernel(n, r0):
    x = np.linspace(-((n-1)//2), (n)//2, n) # Create axis 
    y = np.linspace(-((n-1)//2), (n)//2, n)
    # Create mesh
    xx, yy = np.meshgrid(x, y)
    # Distance from "origin" 
    zz = np.sqrt(xx**2 + yy**2)

    # Apply smoothing
    zz[zz < r0] = r0**2
    zz = 1/zz
    
    return np.roll(np.roll(zz, n//2+1, axis = 0), n//2+1, axis = 1)

def inbound_array_np(xy,n):
    '''Applies periodic conditions but with numpy'''
    xy[xy<-0.5]=xy[xy<-0.5]+n
    xy[xy>=n-0.5]=xy[xy>=n-0.5]-n

@nb.njit
def hist2d_wmass(xy,mat,m):
    '''Compute the density matrix'''
    nx=xy.shape[0] # Number of particles
    for i in range(nx):
        ix=int(xy[i,0]+0.5)
        iy=int(xy[i,1]+0.5)
        if m[i]>0: #we can set m=0 to flag particles
            mat[ix,iy]=mat[ix,iy]+m[i] # Filling up the density matrix
    
def inbound_array_np(xy,n):
    '''Applies periodic conditions but with numpy'''
    xy[xy<-0.5]=xy[xy<-0.5]+n
    xy[xy>=n-0.5]=xy[xy>=n-0.5]-n


def mask_array_np(xy,m,n):
    for i in range(xy.shape[1]):
        m[xy[:,i]<-0.5]=0
        m[xy[:,i]>=n-0.5]=0

class Particles:
    # Class that takes into account all of the particles in our simulation

    def __init__(self, nb_part, nb_grid, soft = 1, periodic = True):
        # Initializing the particles
        

        # Characteristics of particles
        self.x = np.empty([nb_part, 2])     # Positions
        self.f = np.empty([nb_part, 2])     # Forces on particles
        self.v = np.empty([nb_part, 2])     # Velocities
        self.grad = np.empty([nb_part, 2])  # Gradient
        self.m = np.empty(nb_part)          # Mass of particles

        self.kernel = []
        self.kernelft = []
        self.nb_part = nb_part
        self.nb_grid = nb_grid

        if periodic:    # If periodic, just make a normal grid for densities and potentials
            self.rho = np.empty([self.nb_grid, self.nb_grid])     # Creating grid for densities
            self.pot = np.empty([self.nb_grid, self.nb_grid])     # Creating grid for potential
        else:           # If not periodic then make a grid twice as large
            self.rho = np.empty([2*self.nb_grid, 2*self.nb_grid])     # Creating grid for densities
            self.pot = np.empty([2*self.nb_grid, 2*self.nb_grid])     # Creating grid for potential

        self.soft = soft            # Softening factor
        self.periodic = periodic    # Periodic boolean
    
    def get_rho1(self):
        # Function to create the density given the positions of the particles

        # Apply periodic boundaries
        if self.periodic:
            self.x[self.x <= 0] = self.x[self.x <= 0] + self.nb_grid # Lower bound
            self.x[self.x > self.nb_grid] = self.x[self.x > self.nb_grid] - self.nb_grid # Upper bound


        self.rho[:] = 0
        rho_maker(
            (self.x).astype(int),
            self.m,
            self.rho
        )

    def get_rho(self):
        '''Gets density of particles per cell'''
        if self.periodic:
            inbound_array_np(self.x,self.nb_grid)
        else:
            mask_array_np(self.x,self.m,self.nb_grid)
        self.rho[:]=0
        hist2d_wmass(self.x,self.rho,self.m)

    def get_kernel(self):
        # Function to get the kernel
        if self.periodic:   # If space is periodic
            self.kernel = make_kernel(self.nb_grid, self.soft)
        else:               # If space is not periodic
            self.kernel = make_kernel(self.nb_grid*2, self.soft)
        
        # Get FFT of kernel now, so we don't need to do it at every step later
        self.kernelft = np.fft.rfft2(self.kernel)


    def get_pot(self):
        # Function to get the potential using a convolution

        self.get_rho() # Update the density

        n = self.nb_grid

        if not(self.periodic):
            n = n*2 # Make grid larger to avoid "wrap around" feature
        
        # Convolution
        self.pot = np.fft.irfft2(np.fft.rfft2(self.rho) * self.kernelft, [n, n])

    def get_force(self):
        # Function to update the force with the gradient of the potential
        get_grad(self.x, self.pot, self.grad)
        self.f[:] = self.grad
        

    def start_uniform(self, v_scale = 1.0):
        # Create an initial uniformly random distribution of particles
        system.x[:] = np.random.rand(self.nb_part, 2)*self.nb_grid
        system.v[:] = np.random.randn(self.nb_part, 2)*1.0
        system.m[:] = 1
    def ics_corner(self):
        self.x[:] = 1
        self.v[:] = 0
    
    def two_particles(self):
        self.x[0] = np.array([0, 0])
        self.x[1] = np.array([4, 4])

        self.v[0] = np.array([-1, 0])
        self.v[1] = np.array([1, 0])

        self.m[:] = 2

    def update_step(self, dt = 1):
        # Function to make the system take one step forward in time

        # Update position
        self.x[:] = self.x[:] + dt*self.v

        # update potential and density
        self.get_pot()
        self.get_force()

        # Update velocity with force
        self.v[:] = self.v[:] + dt*self.f

n_part = 2
n_grid = 10

system = Particles(
    n_part,
    n_grid
)


#system.start_uniform()
#system.ics_corner()
system.two_particles()

print(system.x)
#system.get_rho()
##system.get_kernel()
#system.get_pot()

print(system.rho)

plt.imshow(system.rho)
plt.show()

assert(0==1)

# Plotting stuff
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
frame = ax.imshow(system.pot)


for i in range(40):
    frame.set_data(system.pot)

    print(system.f)
    system.update_step(1)
    plt.pause(0.001)
