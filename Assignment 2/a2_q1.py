import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
import scipy.constants as cons


def ef_ring(z, Q = 1, R = 1):
    '''Function for electric field of a ring along central axis.'''
    return -Q*z/(4*np.pi*cons.epsilon_0*(z**2+R**2)**(3/2))

v = np.linspace(-4, 4, 1001)
y_data = ef_ring(v)


plt.plot(v, y_data)
plt.title("Electric field along central axis of ring")
plt.xlabel("Distance from center")
plt.ylabel("Electric Field")
plt.savefig("figs/q1_ring_ef.jpg")
plt.show()




