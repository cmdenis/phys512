import numpy as np
import numba as nb
import time
import imageio
from matplotlib import pyplot as plt
from scipy import fft
plt.ion()


def inbound_array_np(xy,n):
    '''Applies periodic conditions with numpy arrays'''
    xy[xy<-0.5] = xy[xy <- 0.5] + n
    xy[xy>=n-0.5] = xy[xy >= n-0.5]-n

def mask_array_np(xy,m,n):
    '''Mask for non-periodic boundary conditions'''
    for i in range(xy.shape[1]):
        m[xy[:,i]<-0.5]=0
        m[xy[:,i]>=n-0.5]=0

@nb.njit(parallel=True)
def get_grad(xy,pot,grad):
    '''Function to compute the gradient of an array'''
    n=pot.shape[0]
    for i in nb.prange(xy.shape[0]):
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
        #potential is f00(1-fx)(1-fy)+f01(1-fx)(fy)+f10(fx)(1-fy)+f11(fx)(fy)
        #grad_x is -f00(1-fy)-f01(fy)+f10(1-fy)+f11(fy)
        #grad_y is -f00(1-fx)+f01(1-fx)-f10(fx)+f11(fx)
        #grad[i,0]=pot[ix0,iy0]*(fy-1)-pot[ix0,iy1]*fy+pot[ix1,iy0]*(1-fy)+pot[ix1,iy1]*ft
        grad[i,0]=(pot[ix1,iy1]-pot[ix0,iy1])*fy+(pot[ix1,iy0]-pot[ix0,iy0])*(1-fy)
        grad[i,1]=(pot[ix1,iy1]-pot[ix1,iy0])*fx+(pot[ix0,iy1]-pot[ix0,iy0])*(1-fx)


@nb.njit
def hist2d_wmass(xy,mat,m):
    nx=xy.shape[0]
    for i in range(nx):
        ix = int(xy[i,0]+0.5)
        iy = int(xy[i,1]+0.5)
        if m[i] > 0: #we can set m=0 to flag particles
            mat[ix,iy]=mat[ix,iy]+m[i]


class particles:
    def __init__(self,npart=10000,n=1000,soft=1,periodic=True):
        self.x = np.empty([npart,2])
        self.f = np.empty([npart,2])
        self.v = np.empty([npart,2])
        self.grad = np.empty([npart,2])
        self.m = np.empty(npart)
        self.kernel = []
        self.kernelft = []
        self.npart = npart
        self.ngrid = n
        if periodic:
            self.rho = np.empty([self.ngrid,self.ngrid])
            self.pot = np.empty([self.ngrid,self.ngrid])
        else:
            self.rho = np.empty([2*self.ngrid,2*self.ngrid])
            self.pot = np.empty([2*self.ngrid,2*self.ngrid])

        self.soft = soft
        self.periodic = periodic


    def two_particles(self):

        n = self.ngrid

        # Positioning the particles at the center of the plane separated from each other
        self.x[0] = np.array([n/2, n*0.4])
        self.x[1] = np.array([n/2, n*0.6])
        factor = 0.5
        self.v[0] = np.array([-1, 0])*factor
        self.v[1] = np.array([1, 0])*factor

        self.m[:] = 8

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
            ngrid2 = ngrid*2
            x=np.fft.fftfreq(ngrid2)*ngrid2
            rsqr=np.outer(np.ones(ngrid2),x**2)
            rsqr=rsqr+rsqr.T
            rsqr[rsqr<soft**2]=soft**2
            self.kernel=rsqr**-0.5

        self.kernelft=fft.rfft2(self.kernel)


    def get_rho(self):
        '''Method to get the density in a grid of the particles.'''
        if self.periodic:
            inbound_array_np(self.x,self.ngrid)
        else:
            mask_array_np(self.x,self.m,self.ngrid)
        self.rho[:]=0
        hist2d_wmass(self.x,self.rho,self.m)
    
    def get_pot(self):
        '''Method to get the potential using a convolution of kernel and density'''
        t1=time.time()
        self.get_rho()
        #print('got density: ',time.time()-t1)
        rhoft=fft.rfft2(self.rho)
        #print('got ft 1: ',time.time()-t1)
        n=self.ngrid
        if not(self.periodic):
            n=n*2
        #self.pot=fft.irfft2(rhoft*self.kernelft,[self.ngrid,self.ngrid])
        self.pot=fft.irfft2(rhoft*self.kernelft,[n,n])
        #print('got ft 2: ',time.time()-t1)
    
    def get_forces(self):
        '''Get the force on every particle using the gradient of the potential.'''
        get_grad(self.x,self.pot,self.grad)
        self.f[:]=self.grad
    
    def take_step(self,dt=1):
        self.x[:]=self.x[:]+dt*self.v
        self.get_pot()
        self.get_forces()
        self.v[:]=self.v[:]+self.f*dt


# Initializing our system with 2 particles
parts=particles(npart=2,n = 50, soft=2,periodic=True)
parts.two_particles()

# Get the kernel
parts.get_kernel()

# Copy 
xy=parts.x.copy()
parts.get_pot()
rho=parts.rho.copy()
pot=parts.pot.copy()

osamp=40    


fig = plt.figure()
ax = fig.add_subplot(111)
crap=ax.imshow(parts.rho[:parts.ngrid,:parts.ngrid]**0.5)

for i in range(500):
    #t1=time.time()
    for j in range(osamp):
        parts.take_step(dt=0.02)
    #t2=time.time()
    #kin=np.sum(parts.v**2)
    #pot=np.sum(parts.rho*parts.pot)
    #print(t2-t1,kin,pot,kin-0.5*pot)
    #assert(1==0)
    #plt.clf()
    #plt.imshow(parts.rho**0.5)#,vmin=0.9,vmax=1.1)
    #plt.colorbar()


    crap.set_data(parts.rho)
    plt.savefig(f'figs/two_particles/{i:003}', dpi = 50)
    plt.pause(0.001)
