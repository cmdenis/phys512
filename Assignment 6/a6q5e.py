import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob


def smooth_vector(vec,sig):
    n=len(vec)
    x=np.arange(n)
    x[n//2:]=x[n//2:]-n
    kernel=np.exp(-0.5*x**2/sig**2) #make a Gaussian kernel
    kernel=kernel/kernel.sum()
    vecft=np.fft.rfft(vec)
    kernelft=np.fft.rfft(kernel)
    vec_smooth=np.fft.irfft(vecft*kernelft) #convolve the data with the kernel
    return vec_smooth

def read_template(filename):
    dataFile=h5py.File(filename,'r')
    template=dataFile['template']
    tp=template[0]
    tx=template[1]
    return tp,tx
def read_file(filename):
    dataFile=h5py.File(filename,'r')
    dqInfo = dataFile['quality']['simple']
    qmask=dqInfo['DQmask'][...]

    meta=dataFile['meta']
    gpsStart=meta['GPSstart'][()]
    utc=meta['UTCstart'][()]
    duration=meta['Duration'][()]
    strain=dataFile['strain']['Strain'][()]
    dt=(1.0*duration)/len(strain)

    dataFile.close()
    return strain,dt,utc

fname = 'LIGO/L-L1_LOSC_4_V2-1135136334-32.hdf5'
tmp_name = 'LIGO/GW150914_4_template.hdf5'

# Loading data
print('reading file ',fname)
strain,dt,utc=read_file(fname)

# Loading template
tp,tx=read_template(tmp_name)

# Creating window function
x=np.linspace(-np.pi/2,np.pi/2,len(strain))
win=np.cos(x)**1


# FT of data
noise_ft=np.fft.rfft(win*strain)
tobs=dt*len(noise_ft)
dnu=1/tobs
nu=np.arange(len(noise_ft))*dnu
nu[0]=0.5*nu[1]

template_ft=np.fft.rfft(tp*win)

noise_smooth=smooth_vector(np.abs(noise_ft)**2, 20)


# Plotting the spectrum for data and template
plt.clf()
plt.loglog(nu[:-1], np.abs(np.abs(noise_ft[:-1])**2/noise_smooth)**2, label = "Whitened Data")
plt.loglog(nu[:-1], np.abs(template_ft[:-1]*np.abs(noise_ft[:-1])**2/noise_smooth)**2, label = "Convolution in Fourier Space")
plt.loglog(nu, np.abs(template_ft)**2, label = "Template")
plt.xlabel("Frequency")
plt.ylabel("Amplitude")
plt.legend()
plt.savefig("figs/a6q5_tmp_and_data.jpg")
plt.show()
