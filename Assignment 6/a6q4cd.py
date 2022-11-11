import numpy as np
import matplotlib.pyplot as plt

N = 100
freq = 10.2
x = np.linspace(0, 2*np.pi, N)
y_analytic = np.sin(x*freq)
window = 0.5 - np.cos(x)/2
w_sine = window*y_analytic

# Taking fft
y_fft = np.fft.rfft(y_analytic)
ft_w_sine = np.fft.rfft(w_sine)
axis = np.arange(len(y_fft))
plt.clf()
plt.plot(np.abs(y_fft), label = "Numerical")

ana = np.zeros(10001)
axis = np.arange(10001)/1000*N/2
ana[axis==10.2]=N/2

# Plotting the analytic and numerical result

plt.plot(axis, ana, label = "Analytical")
plt.legend()
plt.xlabel("Frequency")
plt.ylabel("Amplitude")
plt.savefig("figs/a6q4d_ana_vs_num.jpg")
plt.show()


# Plotting the windowed and unwindowed data
plt.clf()
plt.plot(np.abs(y_fft), label = "Unwindowed")
plt.plot(np.abs(ft_w_sine), label = "Windowed")
plt.legend()
plt.xlabel("Frequency")
plt.ylabel("Amplitude")
plt.savefig("figs/a6q4d_ft_window_comp.jpg")
plt.show()


