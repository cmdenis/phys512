import numpy as np 
import matplotlib.pyplot as plt


def convo_shift(x, shift):
    dirac = np.zeros(len(x))

    dirac[int(len(x)/2)+shift] = 1

    return np.convolve(x, dirac, 'same')

lar = 1000
x = np.arange(lar)

# Making Gaussian
y1 = np.zeros(len(x))
y1[int(len(x)/2)-200] = 1

# Making initial Gaussian
y2 = np.exp(-(x-lar/2)**2/(lar**2/64))


# Applying shifting function through convolution
y3 = convo_shift(y2, -200)


plt.plot(x, y1, label = "Dirac")
plt.plot(x, y2, label = "Gaussian")
plt.plot(x, y3, label = "Convolution")
plt.legend()
plt.savefig("figs/a6q1_convo.jpg")
plt.show()
