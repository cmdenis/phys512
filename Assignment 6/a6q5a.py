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


# Creating windowing function
def windowed(vec):
    x=np.linspace(-np.pi/2,np.pi/2,len(vec))
    win=np.cos(x)
    return vec*win

def getnobel(fname, tmp_name = 'LIGO/GW150914_4_template.hdf5'):

    # Loading data
    print('reading file ',fname)
    strain,dt,utc=read_file(fname)

    template_name= tmp_name
    tp,tx=read_template(template_name)

    x=np.linspace(-np.pi/2,np.pi/2,len(strain))
    win=np.cos(x)

    noise_ft=np.fft.fft(win*strain)
    noise_smooth=smooth_vector(np.abs(noise_ft)**2,10)
    noise_smooth=noise_smooth[:len(noise_ft)//2+1] #will give us same length




    tobs=dt*len(strain)
    dnu=1/tobs
    nu=np.arange(len(noise_smooth))*dnu
    nu[0]=0.5*nu[1]

    Ninv=1/noise_smooth
    Ninv[nu>1500]=0
    Ninv[nu<20]=0

    template_ft=np.fft.rfft(tp*win)
    template_filt=template_ft*Ninv
    data_ft=np.fft.rfft(strain*win)
    rhs=np.fft.irfft(data_ft*np.conj(template_filt))

    plt.plot(rhs)
    plt.show()
    return rhs

plt.plot(getnobel('LIGO/H-H1_LOSC_4_V2-1126259446-32.hdf5'))
plt.show()


# First we load the raw data from Hanford
fname = 'LIGO/H-H1_LOSC_4_V2-1126259446-32.hdf5'
print('reading file ',fname)
strain_h, dt_h, utc_h = read_file(fname)
# FFT of data for Hanford
fft_strain_h = np.fft.rfft(windowed(strain_h))



# We load the data from Livingston
fname = 'LIGO/L-L1_LOSC_4_V1-1167559920-32.hdf5'
print('reading file ',fname)
strain_l, dt_l, utc_l = read_file(fname)
# FFT of data for Livingston
fft_strain_l = np.fft.rfft(strain_l)
# Smoothed data
l_smooth = smooth_vector(np.abs(fft_strain_l)**2,10)

# Defining noise
sig = 0.05



# Loading the template
template_name = 'LIGO/GW150914_4_template.hdf5'
tp, tx = read_template(template_name)
# FFT of template
fft_tp = np.fft.rfft(windowed(tp))
fft_tp = fft_tp*Ninv



# Finding right hand side
rhs = np.fft.irfft(np.conj(fft_tp)*fft_strain_l)/sig**2

# Matched filter noise
mf = rhs/lhs


print('median matched filter noise: ',np.median(np.abs(mf)))
print('median deconvolution noise: ',np.median(np.abs(deconv)))
print('matched filter expected noise: ',1/np.sqrt(lhs))

plt.plot(rhs)
plt.show()

assert(0==1)


# Windowing data
noise_ft=np.fft.rfft(windowed(strain_h))
a_noise_ft = np.abs(noise_ft)


plt.loglog(a_noise_ft, label = "FT of data with windowing")
plt.loglog(smooth_vector(a_noise_ft, 5), label = "FT with smoothing")
plt.xlim([500, 10**5])
plt.legend()
plt.show()

noise_smooth=smooth_vector(np.abs(noise_ft)**2, 5)
noise_smooth=noise_smooth[:len(noise_ft)//2+1] #will give us same length




tobs=dt_l*len(strain_l)
dnu=1/tobs
nu=np.arange(len(noise_smooth))*dnu
nu[0]=0.5*nu[1]

Ninv=1/noise_smooth
Ninv[nu>1500]=0
Ninv[nu<20]=0

template_ft=np.fft.rfft(tp*win)
template_filt=template_ft*Ninv

data_ft=np.fft.rfft(strain_l*win)
rhs=np.fft.irfft(data_ft*np.conj(template_filt))





plt.clf()
plt.loglog(np.abs(rhs))
plt.show()