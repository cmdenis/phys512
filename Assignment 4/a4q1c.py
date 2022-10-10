import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt

# We start by loading the data

stuff = np.load('mcmc/sidebands.npz')
t = stuff['time']
d = stuff['signal']



# We use Newton's method to find the best fit for the data

# We define our fit function
def lorentz_fit(p, t):
    # Evaluating the fit at every point in time
    y = p[0]/(1 + ((t - p[1])/p[2])**2)

    # Initializing the gradient
    grad = np.zeros([t.size, p.size])

    # Finding the derivatives to build the gradient
    grad[:, 0] = 1.0/(1 + ((t - p[1])/p[2])**2)
    grad[:, 1] = 2*p[0]*p[2]**2*(t-p[1])/((t - p[1])**2 + p[2]**2)**2
    grad[:, 2] = 2*p[0]*p[2]*(t-p[1])**2/((t - p[1])**2 + p[2]**2)**2

    return y, grad



times = np.arange(min(t), max(t), t[1]-t[0])

# Initial parameter guess
p = np.array([1, 0.00018, 0.00005])

# Show guess
plt.plot(t, lorentz_fit(p, t)[0], label = "Guess Fit")
plt.plot(t, d, label = "Data")
plt.title("Sideband Raw Signal")



# We do the procedure 5 times to try it out
for j in range(5):
    pred, grad = lorentz_fit(p, t)

    # Residuals
    r = d - pred


    lhs = grad.T@grad
    rhs = grad.T@r 

    dp = np.linalg.inv(lhs)@rhs

    p = p + dp

    print(p, dp)

plt.plot(t, lorentz_fit(p, t)[0], label = "Best Fit")
plt.legend()
plt.savefig("figs/a4q1a_newton_method.jpg")
plt.show()
plt.clf()
    




