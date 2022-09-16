import numpy as np
import matplotlib.pyplot as plt
import scipy as sci

# Function to take the derivative at some point with a given step size
def deriv(func, x0, delta):
    return 2/(3*delta) * ( func(x0 + delta) - func(x0 - delta) - 1/8*(func(x0 + 2*delta) - func(x0 - 2*delta))     )







######## exp(x) Plot

# Log linspace to scan through the different possible stepsizes
dxs = 10**np.linspace(-0.1, -8, 1000)

# Plotting the error as a function of the step size
test_point = 0

plt.loglog(dxs, abs(np.exp(test_point) - deriv(np.exp, test_point, dxs)), label = 'Numerical Error')
plt.plot(10**-3.2, abs(np.exp(test_point) - deriv(np.exp, test_point, 10**-3.2)), 'ro', label = 'Minimum Error Estimation')
plt.title("Error as a function of step size, at x = 0")
plt.xlabel("Step Size")
plt.ylabel("Error")
plt.legend()
plt.savefig("figs/q1_error_plot1.jpg")
plt.show()







####### exp(0.01 Plot)

#Creating the desired function:
def mod_exp(x):
    return np.exp(0.01*x)

dxs = 10**np.linspace(-0, -5, 1000)

# Plotting the error as a function of the step size
plt.loglog(dxs, abs(0.01 - deriv(mod_exp, 0, dxs)), label = 'Numerical Error')
plt.plot(10**-3.2/0.01, abs(0.01 - deriv(mod_exp, 0, 10**-3.2/0.01)), 'ro', label = 'Minimum Error Estimation')
plt.title("Error as a function of step size, at x = 0")
plt.xlabel("Step Size")
plt.ylabel("Error")
plt.legend()
plt.savefig("figs/q1_error_plot2.jpg")
plt.show()


