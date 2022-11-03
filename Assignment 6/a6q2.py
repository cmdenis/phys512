import numpy as np 
import matplotlib.pyplot as plt


# Defining the correlation function
def cor_func(a, b):
    return np.fft.ifft(np.fft.fft(a) * np.conjugate(np.fft.fft(b)))


# Part a)

# Generating our arrays
x = np.linspace(-5, 5, 1001)
y1 = np.exp(-x**2/2)/(np.sqrt(2*np.pi))
y2 = cor_func(y1, y1)
y3 = np.convolve(y1, y1, 'same')

 # plotting
plt.plot(x, y1, label = "Gaussian")
plt.plot(x, y2, label = "Correlation (Fourier space)")
plt.plot(x, y3, label = "Correlation (convolution)")
plt.legend()
plt.savefig("figs/a6q2_a_gauss_cor.jpg")
plt.show()


# Part b)

# Generating our arrays
x = np.linspace(-5, 5, 1001)
y1 = np.exp(-(x-1)**2/2)/(np.sqrt(2*np.pi))
y2 = cor_func(y1, y1)
y3 = np.convolve(y1, y1, 'same')

 # plotting
plt.plot(x, y1, label = "Gaussian")
plt.plot(x, y2, label = "Correlation (Fourier space)")
plt.plot(x, y3, label = "Correlation (convolution)")
plt.legend()
plt.savefig("figs/a6q2_a_gauss_cor.jpg")
plt.show()



