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

# Preparing arrays to be used
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
plt.savefig("figs/a3q1_pt1_comp.jpg")
plt.show()


# Plotting the residuals
plt.title("Residuals")
plt.plot(xs, ys-y_true, label = "RK4")
plt.plot(xs, y_true-y_true, label = "c0 * exp(arctan(x))")
plt.legend()
plt.savefig("figs/a3q1_pt1_resi.jpg")
plt.show()





# Part 2) Using the half step size integrator

def rk4_stepd(fun, x, y, h):
    first_step = rk4_step(fun, x, y, h/2)
    second_step = rk4_step(fun, x + h/2, first_step, h/2)
    return (16*second_step - rk4_step(fun, x, y, h) )/15


# Running the rk4_stepd integrator
ysd = np.array([1])
for x in xs[:-1]:
    ysd = np.append(ysd, rk4_stepd(fun, x, ysd[-1], step_size))


# Comparing the functions

# Plotting the function
plt.plot(xs, ys, label = "rk4_step")
plt.plot(xs, ysd, label = "rk4_stepd")
plt.plot(xs, y_true, label = "c0 * exp(arctan(x))")
plt.legend()
plt.savefig("figs/a3q1_pt2_comp.jpg")
plt.show()


# Plotting the residuals
plt.title("Residuals")
#plt.plot(xs, ys-y_true, label = "rk4_step")
plt.plot(xs, ysd-y_true, label = "rk4_stepd")
plt.plot(xs, y_true-y_true, label = "c0 * exp(arctan(x))")
plt.legend()
plt.savefig("figs/a3q1_pt2_resi.jpg")
plt.show()
