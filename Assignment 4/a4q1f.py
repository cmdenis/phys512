import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt

# We start by loading the data

stuff = np.load('mcmc/sidebands.npz')
t = stuff['time']
d = stuff['signal']

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


times = np.arange(min(t), max(t), t[1]-t[0])

# Initial parameter guess
p = np.array([1.4, 0.0002, 0.00002, 0.2, 0.2, 0.00005])




# We do the procedure 5 times to try it out
for j in range(25):
    # Calculating predicted fit and gradient
    pred = three_lorentz_fit(p, t)
    grad = grad_f(three_lorentz_fit, p, t)
    # Residuals
    r = d - pred
    # Fragmenting matrix operations
    lhs = grad.T@grad
    rhs = grad.T@r 
    # Calculate the new best parameter
    dp = np.linalg.inv(lhs)@rhs
    p = p + dp
    #print(p, dp)

'''Error finding part'''

# Finding the 1) predicted data 2) the grad of the function
pred = three_lorentz_fit(p, t)
grad = grad_f(three_lorentz_fit, p, t)
#print("The parameters are:", p)

# Finding the noise in our data
err = np.mean(np.abs(pred - d))


# Finding the error on the parameters
lhs = grad.T@grad
cov_mat = np.linalg.inv(lhs)*err
p_err = np.sqrt(np.diagonal(cov_mat))

#print("And the error on them are:", p_err)





'''Varying the parameters and looking at the chi^2 part'''

# We define a chi^2 taking function
def chi2(p, x, y):
    pred = three_lorentz_fit(p, x)
    error = np.mean(np.abs(pred - y))
    return np.sum((pred - y)**2/error)

print(cov_mat.shape[1])

alt_chi = np.zeros(len(p))
plt.plot(t, d, label = "Data")
for i in range(len(p)):
    pert = np.zeros(len(p)) 
    pert[i] = np.sqrt(cov_mat[i, i])
    alt_chi[i] = chi2(p + pert, t, d) - chi2(p, t, d)
    plt.plot(t, three_lorentz_fit(p + pert, t), label = "Param perturbation "+str(i+1))



print([np.sqrt(cov_mat[i, i]) for i in range(6)])
print(alt_chi)










'''Plotting part'''


plt.title("Sideband Residuals")
plt.xlabel("Time (t)")
plt.ylabel("Amplitude")
plt.plot(t, three_lorentz_fit(p, t), label = "Best Fit")
plt.legend()
plt.savefig("figs/a4q1e_residuals.jpg")
plt.show()
plt.clf()
    




