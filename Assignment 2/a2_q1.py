from ast import arg
import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
import scipy.constants as cons
from scipy import integrate

# We define our own function for integration:
def integrate_adaptive(fun, a, b, tol, extra = None):
    print('calling adaptive function from ', a, b)
    sub_div = 5
    if extra == None:
        x = np.linspace(a, b, sub_div)
        dx = x[1] - x[0]
        y = fun(x)
    else:
        x = np.linspace(a, b, sub_div)
        dx = x[1] - x[0]
        y = np.array([extra[0], fun(x[1]), extra[1], fun(x[3]), extra[2]])

    #do the 3-point integral
    i1 = (y[0]+4*y[2]+y[4])/3*(2*dx)
    i2 = (y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/3*dx
    myerr = np.abs(i1-i2)
    if myerr < tol:
        return i2
    else:
        mid = (a+b)/2
        int1 = integrate_adaptive(fun, a, mid, tol/2, extra = [y[0], y[1], y[2]])
        int2 = integrate_adaptive(fun, mid, b, tol/2, extra = [y[2], y[3], y[4]])
        return int1+int2

R = 1
sigma = 0.1*(4 * np.pi * R**2)**-1

def to_int(theta, z):
    coeff = 1#2*np.pi*R**2*sigma/(4*np.pi*cons.epsilon_0)
    return -coeff*(z-R*np.cos(theta))*np.sin(theta)/(R**2 + z**2 - 2*R*z*np.cos(theta))**(3/2)

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

