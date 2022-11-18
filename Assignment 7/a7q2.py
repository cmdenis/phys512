import numpy as np
import matplotlib.pyplot as plt


# Comparison between the different bounding distributions
def power_law(x):
    return 1/(x+1)#/np.log(2)
def gauss(x):
    return np.exp(-x**2/2)#/(2 - 2/np.sqrt(np.e))
def cauchy(x):
    return 1/(x**2 + 1)#*4/np.pi
def expo(x):
    return np.exp(-x)#*np.e/(np.e - 1)

x = np.linspace(0, 5, 1001)

plt.plot(x, expo(x), label = "Exponential")
plt.plot(x, power_law(x), label = "Power Law")
plt.plot(x, gauss(x), label = "Gaussian")
plt.plot(x, cauchy(x), label = "Cauchy")
plt.plot(x, x*0+1, label = "Uniform")
plt.title("Comparison between different distributions")
plt.legend()
plt.savefig("figs/a7q2_comp_dist.jpg")
#plt.show()
plt.clf()


# First we define a function to get a cauchy distribution and make sure it works

# Using a recursive function
def cauchy_rand():
    # Generating 2 random numbers
    r = np.random.rand(2)*4
    r[0] = r[0]
    # Check if they are accepted, else try again recursively
    if r[1] < cauchy(r[0]):
        return r[0]
    else:
        return cauchy_rand()


cauchy_samples = np.array([cauchy_rand() for i in range(100000)])
bins = np.linspace(0, 4, 50)
plt.hist(cauchy_samples, bins, density = True, label = "Samples From Cauchy")
plt.plot(bins, cauchy(bins)/np.arctan(4), label = "Cauchy Distribution")
plt.plot(bins, expo(bins)*(1 - np.cosh(4) + np.sinh(4)), label = "Exponential Distribution")
plt.legend()
plt.savefig("figs/a7q2_cauchy.jpg")
plt.show()
plt.clf()




def rand_exp(n):
    # Returns distribution with exponential distribution

    # Compute two random numbers
    r1 = np.random.rand(n)*4
    r2 = np.array([cauchy_rand() for i in range(n)])


    accept = r1 < expo(r2)/cauchy(r2)

    return r2[accept]

n = 1000000
bins = np.linspace(0, 4, 50)
rand_samples = rand_exp(n)
print("The ratio of acceptance is:", len(rand_samples)/n)

plt.clf()
plt.hist(rand_samples, bins, density = True, label = "Samples")
plt.plot(bins, expo(bins)/(1 - np.cosh(4) + np.sinh(4)), label = "Exponential Distribution")
plt.legend()
plt.savefig("figs/a7q2_expo_dist.jpg")
plt.show()
plt.clf()


