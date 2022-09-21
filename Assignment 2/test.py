from ast import Lambda, arg
import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
import scipy.constants as cons
from scipy import integrate

R = 1
sigma = 0.1*(4 * np.pi * R**2)**-1

def to_int(theta, z):
    coeff = 2*np.pi*R**2*sigma/(4*np.pi*cons.epsilon_0)
    return coeff*(z-R*np.cos(theta))*np.sin(theta)/(R**2 + z**2 - 2*R*z*np.cos(theta))**(3/2)

#print(integrate.quad(to_int, -np.pi/2, np.pi/2, args=(5,)))

def efield(z):
    temp_func = lambda theta : to_int(theta, z) # Anonymous function

    return integrate.quad(temp_func, 0, np.pi)

print(0.1/(4*np.pi*cons.epsilon_0*1))


v = np.linspace(-4, 4, 1001)
y_data = [efield(i)[0] for i in v]
#y_coul = 0.1/(4*np.pi*cons.epsilon_0*v**2)

plt.clf()
plt.plot(v, y_data, label = "")
#plt.plot(v, y_coul)
plt.show()