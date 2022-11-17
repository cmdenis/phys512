import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-4, 4, 1001)


def gaussian(x, mean, sigma):
    denominator = np.sqrt(2*np.pi)*sigma
    return np.exp(-((x-mean)/sigma)**2/2)/denominator

sigma = 1
mean = 0
gauss = gaussian(x, mean, sigma)

N = 100001
# Array of random number
rand_n = np.zeros(N)



# Using a recursive function
def gaussian_rand(mean, sigma):
    # Generating 2 random numbers
    r = np.random.rand(2)
    r[0] = r[0]*20 - 10
    # Check if they are accepted, else try again recursively
    if r[1] < gaussian(r[0], mean, sigma):
        return r[0]
    else:
        return gaussian_rand(mean, sigma)




# Using vectorizing and numpy arrays
def rand_np(n, range, pdf):
    dist = range[1] - range[0] # Distance between highest and lowest value in desired range
    r1 = np.random.rand(n)*dist + range[0]
    r2 = np.random.rand(n)

    accept = r2 < pdf(r1)

    return r1[accept]

def f(x):
    sigma = 1
    mean = 0
    denominator = np.sqrt(2*np.pi)*sigma
    return np.exp(-((x-mean)/sigma)**2/2)/denominator

gaussian_data = rand_np(1000000, [x[0], x[-1]], f)

#print(len(gaussian_data))

bins = np.linspace(x[0], x[-1], 50)
plt.hist(gaussian_data, bins, density = True)
plt.plot(x, f(x))
#plt.show()
plt.clf()





# Using vectorizing and numpy arrays
def PDF(x):
    return

def rand_exp(n, range):

    dist = range[1] - range[0] # Distance between highest and lowest value in desired range
    r1 = np.random.rand(n)*dist + range[0]
    r2 = np.random.rand(n)

    accept = r2 < pdf(r1)

    return r1[accept]

bins = np.linspace(-10, 10, 50)
plt.hist(gaussian_data, bins, density = True)
plt.plot(x, f(x))
#plt.show()
plt.clf()


# Making a power law shaped distribution
def power_law(x):
    return 1/(x+1)/np.log(2)

# Using a recursive function
def pl_rand():
    # Generating 2 random numbers
    r = np.random.rand(2)
    r[0] = r[0]
    # Check if they are accepted, else try again recursively
    if r[1] < power_law(r[0])*np.log(2):
        return r[0]
    else:
        return pl_rand()

x = np.linspace(0, 1, 1001)
power_law_samples = np.array([pl_rand() for i in range(100000)])
bins = np.linspace(0, 1, 50)
plt.hist(power_law_samples, bins, density = True)
plt.plot(x, power_law(x))
plt.show()

