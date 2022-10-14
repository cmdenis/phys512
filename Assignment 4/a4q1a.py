import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt



# We start by loading the data
stuff = np.load('mcmc/sidebands.npz')
t = stuff['time']
d = stuff['signal']

# Define a chi^2 function
def chi2(p, x, y):
    pred, grad = lorentz_fit(p, x)
    error = np.mean(np.abs(pred - y))
    return np.sum((pred - y)**2/error)

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


# Initial parameter guess
p = np.array([1, 0.00018, 0.00005])

# Show guess
plt.plot(t, lorentz_fit(p, t)[0], label = "Guess Fit")
plt.plot(t, d, label = "Data")
plt.title("Sideband Signal with Fit")

# Starting the initial loop conditions
looper = True
old_chi = chi2(p, t, d)

# We do the procedure 10 times to try it out
while looper:
    # Find useful prediction and grad matrix
    pred, grad = lorentz_fit(p, t)

    # Residuals
    r = d - pred

    # Fragmenting linear algera stuff
    lhs = grad.T@grad
    rhs = grad.T@r 

    # Applying a step
    dp = np.linalg.inv(lhs)@rhs
    p = p + dp

    new_chi = chi2(p, t, d)

    # Check if if Chi^2 is minimized enough
    if (old_chi - new_chi) < 0.01:
        looper = False
    else:
        old_chi = new_chi

    #print("\nParameters:", p, "\nDisplacement", dp)


# Plot the results
plt.plot(t, lorentz_fit(p, t)[0], label = "Best Fit")
plt.legend()
plt.xlabel("Time (t)")
plt.ylabel("Amplitude")
plt.savefig("figs/a4q1a_newton_method.jpg")
plt.show()
plt.clf()
    
# Print the results
print("The best-fit parameters are:", p)
