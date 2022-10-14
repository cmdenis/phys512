import numpy as np 
from tqdm import tqdm
from scipy import integrate
from matplotlib import pyplot as plt
import corner

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

# Show guess

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

# Error finding part

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




'''Now the MCMC part'''

step_n = 10000                                 # Number of steps to take

p_init = p.copy()*1                          # Initial Parameters

step_size = np.sqrt(np.diagonal(cov_mat))   # Size of steps

chain = np.zeros([step_n, len(p_init)+1])   # Initializing the chain
chain[0, 0:-1] = p_init

cur_pos = chain[0, 0:-1]

cur_chi = chi2(cur_pos, t, d)

chain[0, -1] = cur_chi


test = step_size*np.random.randn(len(step_size))
#print(test)
#print(chi2(chain[0, 0:-1]+test, t, d))

#assert(0==1)
for i in tqdm(range(1, step_n)):

    # Finding new position
    new_pos = chain[i-1, 0:-1] + step_size*np.random.randn(len(step_size))

    # Finding new chi square
    new_chi = chi2(new_pos, t, d)

    #print("The new Chi^2 is:", new_chi)
    #print("The old Chi^2 is:", cur_chi)

    if new_chi < cur_chi:
        #print("Then we accept the change.")
        accept = True
    else:
        delt = new_chi - cur_chi
        #print("Then we don't immediately accept the change. Their difference is", delt)
        
        prob = np.exp(-0.5*delt)
        if np.random.rand() < prob:
            accept = True
        else:
            accept = False
    if accept:
        cur_chi = new_chi
        cur_pos = new_pos

    chain[i, 0:-1] = cur_pos
    chain[i, -1] = cur_chi


# Printing values
best_param = np.mean(chain[4000:-1, :], axis = 0)
param_std = np.std(chain[4000:-1, :], axis = 0)

print("The best fit parameters are:", best_param[0:-1])
print("Their error is:", param_std[0:-1])
print("With Chi^2:", chi2(best_param[0:-1], t, d))

plt.plot(chain[:, 0:-1], label = ["a", "t_0", "w", "b", "c", "dt"])
plt.title("Parameters Trajectory in MCMC")
plt.legend()
plt.xlabel("MCMC Steps")
plt.ylabel("Parameter Value")
plt.savefig("figs/a4q1g_mcmc_parameters.jpg")
plt.show()
plt.clf()

corner.corner(
    chain[:, 0:-1],
    labels=[
        r"$a$",
        r"$t_0$",
        r"$w$",
        r"$b$",
        r"$c$",
        r"$dt$",
    ]
    )
plt.savefig("figs/a4_mcmc_cornerplot.jpg")

plt.show()

          

    




