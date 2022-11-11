import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
import json



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
    #gpsStart=meta['GPSstart'].value
    gpsStart=meta['GPSstart'][()]
    #print meta.keys()
    #utc=meta['UTCstart'].value
    utc=meta['UTCstart'][()]
    #duration=meta['Duration'].value
    duration=meta['Duration'][()]
    #strain=dataFile['strain']['Strain'].value
    strain=dataFile['strain']['Strain'][()]
    dt=(1.0*duration)/len(strain)

    dataFile.close()
    return strain,dt,utc



#fnames=glob.glob("[HL]-*.hdf5")
#fname=fnames[0]
fname='LIGO/H-H1_LOSC_4_V2-1128678884-32.hdf5'
print('reading file ',fname)
strain,dt,utc=read_file(fname)

#th,tl=read_template('GW150914_4_template.hdf5')
template_name='LIGO/LVT151012_4_template.hdf5'
tp,tx=read_template(template_name)

x=np.linspace(-np.pi/2,np.pi/2,len(strain))
win=np.cos(x)

noise_ft=np.fft.fft(win*strain)
print(len(noise_ft))

plt.loglog(np.abs(noise_ft)**2)


noise_smooth=smooth_vector(np.abs(noise_ft)**2,20)
#noise_smooth=noise_smooth[:len(noise_ft)//2+1] #will give us same length
w_data = noise_ft/noise_smooth
plt.loglog(noise_smooth)

plt.loglog(np.abs(noise_ft/noise_smooth)**2)

plt.show()

# Match filtering
template_ft=np.fft.fft(tp*win)
rhs=np.fft.irfft(w_data*np.conj(template_ft))

plt.plot(rhs)
plt.show()


assert(0==1)
tobs=dt*len(strain)
dnu=1/tobs
nu=np.arange(len(noise_smooth))*dnu
nu[0]=0.5*nu[1]

Ninv=1/noise_smooth
Ninv[nu>1500]=0
Ninv[nu<20]=0
print(len(np.fft.rfft(strain)))

template_ft=np.fft.rfft(tp*win)
template_filt=template_ft*Ninv
data_ft=np.fft.rfft(strain*win)
print(len(data_ft))
rhs=np.fft.irfft(data_ft*np.conj(template_filt))

plt.plot(rhs)
plt.show()