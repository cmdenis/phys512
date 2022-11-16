import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-10, 10, 1001)


def gaussian(x, mean, sigma):
    denominator = np.sqrt(2*np.pi)*sigma
    return np.exp(-((x-mean)/sigma)**2/2)/denominator

sigma = 1
mean = 0
gauss = gaussian(x, mean, sigma)

N = 100001
# Array of random number
rand_n = np.zeros(N)




def gaussian_rand(mean, sigma):
    # Generating 2 random numbers
    r = np.random.rand(2)
    r[0] = r[0]*20 - 10

    # Check if they are accepted, else try again recursively
    if r[1] < gaussian(r[0], mean, sigma):
        return r[0]
    else:
        return gaussian_rand(mean, sigma)

gaussian_data = np.array([gaussian_rand(mean, sigma) for i in range(N)])



bins = np.linspace(-5, 5, 51)
plt.hist(gaussian_data, bins, density = True)
plt.plot(x, gauss)
plt.show()



