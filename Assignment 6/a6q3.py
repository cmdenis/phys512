import numpy as np 
import matplotlib.pyplot as plt

# Defining the convolution function
def convo(a, b):
    # Adding zeros to the end
    nb_zer = len(a)
    zers = np.zeros(nb_zer)
    a_zer = np.concatenate([a, zers])
    b_zer = np.concatenate([b, zers])

    # Taking the rFFT and then converting back into position space
    return np.fft.irfft(np.fft.rfft(a_zer) * np.conjugate(np.fft.rfft(b_zer)))[:nb_zer]


# Simple test to see if this is working
# Generating our arrays
x1 = np.array([1, 1, 4, 1, 1, 1])
x2 = np.array([1, 2, 0, 0, 0, 0])

print(convo(x1, x2))

