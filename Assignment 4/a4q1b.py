import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt

'''We start by loading the data'''

stuff = np.load('mcmc/sidebands.npz')
t = stuff['time']
d = stuff['signal']



'''We use Newton's method to find the best fit for the data'''

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





'''We do the "Newton method" 10 times'''

nb_iter = 10

# Initial parameter guess
p = np.array([1, 0.00018, 0.00005])
plt.plot(t, lorentz_fit(p, t)[0], label = "Guess Fit")

# Loop for each step in Newton's method
for j in range(nb_iter):

    # Computing 1) Predicted points 2) Gradient
    pred, grad = lorentz_fit(p, t)

    # Residuals
    r = d - pred

    # Useful matrices
    lhs = grad.T@grad
    rhs = grad.T@r 

    # Computing the step
    dp = np.linalg.inv(lhs)@rhs

    # Finding a closer parameter
    p = p + dp




'''Error finding part'''

# Finding the 1) predicted data 2) the grad of the function
pred, grad = lorentz_fit(p, t)
print("The parameters are:", p)

# Finding the noise in our data
err = np.mean(np.abs(pred - d))


# Finding the error on the parameters
lhs = grad.T@grad
cov_mat = np.linalg.inv(lhs)*err
p_err = np.sqrt(np.diagonal(cov_mat))

print("And the error on them are:", p_err)

plt.plot(t, d, label = "Data")
plt.plot(t, lorentz_fit(p, t)[0], label = "Best Fit")
plt.legend()
plt.title("Sideband Signal With Fits")
plt.savefig("figs/a4q1a_newton_method.jpg")
plt.show()
plt.clf()