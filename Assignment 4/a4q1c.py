import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt

# We start by loading the data
stuff = np.load('mcmc/sidebands.npz')
t = stuff['time']
d = stuff['signal']

# Define a chi^2 function
def chi2(p, x, y):
    pred = lorentz_fit(p, x)
    error = np.mean(np.abs(pred - y))
    return np.sum((pred - y)**2/error)

# First we create a derivative taking function
def p_deriv(func, p_ind, p, t):
    shift = 1e-8
    # Creating array
    dp = np.zeros(len(p))
    # Check if derivative index is an integer
    if isinstance(p_ind, int):
        dp[p_ind] = shift 
    else:
        raise ValueError("Derivative index must be an integer.")
    return (func(p + dp, t) - func(p - dp, t))/(2*shift)    # Two sided derivative

# We use Newton's method to find the best fit for the data

# We define our fit function
def lorentz_fit(p, t):
    # Lorentz function
    return p[0]/(1 + ((t - p[1])/p[2])**2)

# Define numerical gradient taker
def grad_f(f, p, t):
    # Initializing the gradient
    grad = np.zeros([t.size, p.size])
    # Finding the derivatives to build the gradient
    for param in range(len(p)):
        grad[:, param] = p_deriv(f, param, p, t)
    return grad

times = np.arange(min(t), max(t), t[1]-t[0])

# Initial parameter guess
p = np.array([1, 0.00018, 0.00005])

# Show guess
plt.plot(t, lorentz_fit(p, t), label = "Guess Fit")
plt.plot(t, d, label = "Data")
plt.title("Signal with fit")

# Starting the initial loop conditions
looper = True
old_chi = chi2(p, t, d)

# We do the procedure until chi^2 condition is satisfied
while looper:
    # Calculating predicted fit and gradient
    pred = lorentz_fit(p, t)
    grad = grad_f(lorentz_fit, p, t)
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


plt.plot(t, lorentz_fit(p, t), label = "Best Fit")
plt.legend()
plt.xlabel("Time (t)")
plt.ylabel("Amplitude")
plt.annotate("The parameters are: "+str(p), (0.000, 1.4), fontsize = 6)
plt.savefig("figs/a4q1c_newton_method.jpg")
plt.show()
plt.clf()
    
print("The best-fit parameters are:", p)
