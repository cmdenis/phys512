import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt

# We start by loading the data
stuff = np.load('mcmc/sidebands.npz')
t = stuff['time']
d = stuff['signal']

# Define a chi^2 function
def chi2(p, x, y):
    pred = three_lorentz_fit(p, x)
    error = np.mean(np.abs(pred - y))
    return np.sum((pred - y)**2/error)

# First we create a derivative taking function
def p_deriv(func, p_ind, p, t):
    shift = 1e-6
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
def three_lorentz_fit(p, t):
    # Lorentz function
    return p[0]/(1 + ((t - p[1])/p[2])**2) + p[3]/(1 + ((t - p[1]+p[5])/p[2])**2) + p[4]/(1 + ((t - p[1]-p[5])/p[2])**2)

# Define numerical gradient taker
def grad_f(f, p, t):
    # Initializing the gradient
    grad = np.zeros([t.size, p.size])
    # Finding the derivatives to build the gradient
    for param in range(len(p)):
        grad[:, param] = p_deriv(f, param, p, t)
    return grad

# Initial parameter guess
p = np.array([1.4, 0.00019, 0.00002, 0.15, 0.15, 0.00005])

# Starting the initial loop conditions
looper = True
old_chi = chi2(p, t, d)
trunc = 0

# We do the procedure until chi^2 condition is satisfied
while looper:
    # Calculating predicted fit and gradient
    pred = three_lorentz_fit(p, t)
    grad = grad_f(three_lorentz_fit, p, t)
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


'''Error finding part'''

# Estimating the noise in our data
err = np.std(d[-100:-1])
print("Estimated noise:", err)

# Finding the 1) predicted data 2) the grad of the function
pred = three_lorentz_fit(p, t)
grad = grad_f(three_lorentz_fit, p, t)
print("The parameters are:", p)

# Finding the error on the parameters
lhs = grad.T@grad
cov_mat = np.linalg.inv(lhs)*err
p_err = np.sqrt(np.diagonal(cov_mat))

print("And the error on them are:", p_err)

# Plotting part
plt.plot(t, d - three_lorentz_fit(p, t), label = "Data")
plt.plot(t, 0*three_lorentz_fit(p, t), label = "Best Fit")
plt.title("Newton Method Residuals")
plt.legend()
plt.xlabel("Time (t)")
plt.ylabel("Amplitude")
plt.savefig("figs/a4q1e_newton_method.jpg")
plt.show()
plt.clf()
