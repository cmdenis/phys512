import numpy as np
import matplotlib.pyplot as plt

# Number of samples to try
n = 100000

# Define exponential
def exponential(x):
    return np.exp(-x)

# Generate two uniform distribution with proper scaling
u = np.random.rand(n)
v = np.random.rand(n)*2/np.e

# Generate accept array
accept = u <= np.sqrt(exponential(v/u))

# Select accepted values
values = (v/u)[accept]

print("Ratio of acceptance:", len(values)/len(accept))

# Plotting
bins = np.linspace(0, 4, 50)
plt.hist(values, bins=bins, density=True, label = "Samples")
plt.plot(bins, exponential(bins)/(1 - np.cosh(4) + np.sinh(4)), label = "Exponential PDF")
plt.legend()
plt.savefig("figs/a7q3_hist.jpg")
plt.show()




