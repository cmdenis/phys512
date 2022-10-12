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
    return (func(p + dp, t) - func(p, t))/shift


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
for j in range(5):
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


plt.plot(t, d- three_lorentz_fit(p, t), label = "Data")
plt.title("Sideband Residuals")
plt.xlabel("Time (t)")
plt.ylabel("Amplitude")
plt.plot(t, three_lorentz_fit(p, t)*0, label = "Best Fit")
plt.legend()
plt.savefig("figs/a4q13_residuals.jpg")
#plt.show()
plt.clf()

def chi2(p, x, y):
    pred = three_lorentz_fit(p, x)
    error = np.mean(np.abs(pred - y))
    return np.sum((pred - y)**2/error)




step_n = 10000                                 # Number of steps to take

p_init = p.copy()                           # Initial Parameters

step_size = np.array([0.00000001, 0.000000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001])   # Size of steps

chain = np.zeros([step_n, len(p_init)+1])   # Initializing the chain
chain[0, 0:-1] = p_init

cur_pos = chain[0, 0:-1]

cur_chi = chi2(cur_pos, t, d)

chain[0, -1] = cur_chi


test = step_size*np.random.randn(len(step_size))
#print(test)
#print(chi2(chain[0, 0:-1]+test, t, d))

#assert(0==1)
for i in range(1, step_n):

    # Finding new position
    new_pos = chain[i-1, 0:-1] + step_size*np.random.randn(len(step_size))

    # Finding new chi square
    new_chi = chi2(new_pos, t, d)

    print("The new Chi^2 is:", new_chi)
    print("The old Chi^2 is:", cur_chi)

    if new_chi < cur_chi:
        print("Then we accept the change.")
        accept = True
    else:
        delt = new_chi - cur_chi
        print("Then we don't immediately accept the change. Their difference is", delt)
        
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


#print(chain[:, 0:-1])

plt.plot(chain[:, 0])
plt.show()

          

    




