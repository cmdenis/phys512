import numpy as np 
import matplotlib.pyplot as plt

# Defining the convolution function
def convo(a, b):
    # Adding zeros to the end
    nb_zer = 20
    zers = np.zeros(nb_zer)
    a_zer = np.concatenate([zers, a, zers])
    b_zer = np.concatenate([zers, b, zers])
    return np.fft.ifft(np.fft.fft(a_zer) * np.conjugate(np.fft.fft(b_zer)))[nb_zer-1:-nb_zer-1]

# Generating our arrays
x = np.linspace(-5, 5, 1001)
y1 = np.exp(-(x-1)**2/2)/(np.sqrt(2*np.pi))
y2 = convo(y1, y1)
y3 = np.convolve(y1, y1, 'same')


 # plotting
plt.plot(x, y1, label = "Gaussian")
plt.plot(x, y2, label = "Correlation (Fourier space)")
plt.plot(x, y3, label = "Correlation (Convolution)")
plt.legend()
plt.savefig("figs/a6q3.jpg")
plt.show()
