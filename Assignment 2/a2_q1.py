from ast import arg
import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
import scipy.constants as cons
from scipy import integrate

R = 1 # We define the radius of the ball to be 1
Q = (4 * np.pi * R**2)**-1 # We define the surface charge density


def ef_ring(r, z):
    '''Function for electric field of a ring along central axis.'''
    return -Q*z/(4*np.pi*cons.epsilon_0*(z**2+r**2)**(3/2))

# Since we have a ball, let's assume the test charge is at a distance z away from the ball
# As we change the radius, the rings get closer and smaller from that point that is not moving
# We must find the relationship between the radius of the ring and its "z" coordinate, it is:
def new_z(r):
    return np.sqrt(1-r/R)

def int_rings(theta, z1):

    z = z1 - np.sin(theta) # To account for proximity of rings
    r = np.cos(theta)
    return -Q*z/(4*np.pi*cons.epsilon_0*(z**2+r**2)**(3/2))*(2*np.pi*R)*np.sin(theta)*R



v = np.linspace(-4, 4, 1001)
y_data = ef_ring(1, v)


plt.plot(v, y_data)
plt.title("Electric field along central axis of ring")
plt.xlabel("Distance from center")
plt.ylabel("Electric Field")
plt.savefig("figs/q1_ring_ef.jpg")
#plt.show()



# Integrating all rings from with scipy.quad

def test(x, a):
    return a*np.exp(x)


def quad_shell_ef(z):

    return integrate.quad(ef_ring, 1)


def integrand(a, b): 
    #print(theta, t)
    return np.exp(a)*0 + b

#print(integrate.quad(test, 0, 1), args = (1,))

print(integrate.quad(integrand, 0, 1, args=(3,)))

print(integrate.quad(int_rings, -R, R, args=(5,)))

