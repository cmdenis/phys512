import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
import scipy.constants as cons
from scipy import integrate

# We define our own function for integration (based on Question 2):
def integrate_adaptive(fun, a, b, tol, extra = None):
    # print('calling adaptive function from ', a, b)
    sub_div = 5
    if extra == None:
        x = np.linspace(a, b, sub_div)
        dx = x[1] - x[0]
        y = fun(x)
    else:
        x = np.linspace(a, b, sub_div)
        dx = x[1] - x[0]
        y = np.array([extra[0], fun(x[1]), extra[1], fun(x[3]), extra[2]])

    # Approximate Integral
    i1 = (y[0]+4*y[2]+y[4])/3*(2*dx)
    i2 = (y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/3*dx
    myerr = np.abs(i1-i2)
    if myerr < tol:
        return i2
    else:
        mid = (a+b)/2
        i1 = integrate_adaptive(fun, a, mid, tol/2, extra = [y[0], y[1], y[2]])
        i2 = integrate_adaptive(fun, mid, b, tol/2, extra = [y[2], y[3], y[4]])
        return i1+i2





R = 1 # Radius of sphere
sigma = 0.1*(4 * np.pi * R**2)**-1 # Surface charge density (total is 1 C)

def to_int(theta, z):
    # Function to integrate. Without loss of generality we rescale it by making the coefficient 1.
    coeff = 1 #2*np.pi*R**2*sigma/(4*np.pi*cons.epsilon_0)
    return -coeff*(z-R*np.cos(theta))*np.sin(theta)/(R**2 + z**2 - 2*R*z*np.cos(theta))**(3/2)


# Plotting the function to integrate
x = np.linspace(0, np.pi, 100)
y2 = to_int(x, 1.1)
y3 = to_int(x, 3)
y4 = to_int(x, 4)

plt.plot(x, y2, label = "z = 1.1")
plt.plot(x, y3, label = "z = 3")
plt.plot(x, y4, label = "z = 4")
plt.xlabel("Distance from center of sphere")
plt.legend()
plt.title("Function to integrate")
plt.savefig("figs/a2q1_func_int.jpg")
plt.show()
plt.clf()



# Defining the electric field functions (with the built in integration)

def efield_quad(z):
    '''Function that integrates the above function for a given z. Done with quad.'''
    temp_func = lambda theta : to_int(theta, z) # Anonymous function
    return integrate.quad(temp_func, 0, np.pi) # Integration


def efield_custom(z):
    '''Function that integrates the above function for a given z. Done with custom adaptive integrator.'''
    temp_func = lambda theta : to_int(theta, z) # Anonymous function
    return integrate_adaptive(temp_func, 0.00001, np.pi, 0.001) # Integration



# Plotting the electric field with the two integrators

v = np.linspace(0, 4, 1001) # Linspace
y_data_quad = [efield_quad(i)[0] for i in v]  # y data for quad
y_data_custom = [efield_custom(i) for i in v] # y data for custom integrator
#y_coul = 0.1/(4*np.pi*cons.epsilon_0*v**2)

plt.clf()
plt.plot(v, y_data_quad, label = "integrate.quad")
plt.plot(v, y_data_custom, label = "integrate_adaptive")
plt.legend()
plt.ylabel("Electric Field")
plt.xlabel("Distance from center of sphere")
plt.title("Electric Field as function of the distance from the ball")
plt.savefig("figs/a2q1_efield_ball.jpg")
plt.show()




