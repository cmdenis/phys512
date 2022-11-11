import numpy as np
import matplotlib.pyplot as plt

N = 50
freq = 10.2
x = np.linspace(0, 2*np.pi, N)
y_analytic = np.sin(x*freq)
window = 0.5 - np.cos(x)/2
w_sine = window*y_analytic

# Plotting the fourier transform of the window
ft_win = np.fft.fft(window)
axis = np.arange(len(ft_win))
plt.clf()
plt.plot(axis, np.abs(ft_win))
plt.xlabel("Frequency")
plt.ylabel("Amplitude")
plt.savefig("figs/a6q4e_ft_window.jpg")
plt.show()


# Taking fft
y_fft = np.fft.rfft(y_analytic)
ft_w_sine = np.fft.rfft(w_sine)
ft_win = np.fft.rfft(window)
axis = np.arange(len(ft_win))

# Manually convolving the signal the fourier transform of the sine with a "neighbor coupling"
convo_spec = []
for i, j in enumerate(y_fft):
    convo_spec.append(
        N/2*j - N/4*(y_fft[i-1] + y_fft[np.mod(i+1, len(y_fft))])
        )
convo_spec = np.array(convo_spec)/N


# Plotting the convolved and related spectra
plt.clf()
plt.plot(np.abs(y_fft), label = "Unwindowed")
plt.plot(np.abs(convo_spec), label = "Convolved")
plt.plot(np.abs(ft_w_sine), label = "Windowed")
plt.legend()
plt.xlabel("Frequency")
plt.ylabel("Amplitude")
plt.savefig("figs/a6q4e_ft_window_comp.jpg")
plt.show()


