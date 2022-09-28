import numpy as np
from matplotlib import pyplot as plt 


# Part 1) evaluating the function using rk4

# Defining single step integrator rk4
def rk4_step(fun, x, y, h):
    '''Runge-Kutta 4 step integrator'''
    k1 = h*fun(x, y)
    k2 = h*fun(x+h/2, y+k1/2)
    k3 = h*fun(x+h/2, y+k2/2)
    k4 = h*fun(x+h, y+k3)
    return y + (k1 +2*k2 + 2*k3 +k4)/6


# Defining function to use
fun = lambda x, y : y/(1+x**2)

xs = np.linspace(-20, 20, 200)
step_size = xs[1] - xs[0]
ys = np.array([1])


# Evaluating the function at all the time steps
for x in xs[:-1]:
    ys = np.append(ys, rk4_step(fun, x, ys[-1], step_size))


# Creating array with the "true" values
c0 = np.exp(np.arctan(20))
y_true = c0 * np.exp(np.arctan(xs))

# Plotting the function
plt.plot(xs, ys, label = "RK4")
plt.plot(xs, y_true, label = "c0 * exp(arctan(x))")
plt.legend()
plt.show()


# Plotting the residuals
plt.title("Residuals")
plt.plot(xs, ys-y_true, label = "RK4")
plt.plot(xs, y_true-y_true, label = "c0 * exp(arctan(x))")
plt.legend()
plt.show()


