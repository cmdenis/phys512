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



def getnobel(fname, tmp_name = 'LIGO/GW150914_4_template.hdf5'):

    # Loading data
    print('reading file ',fname)
    strain,dt,utc=read_file(fname)

    # Loading template
    tp,tx=read_template(tmp_name)

    # Creating window function
    x=np.linspace(-np.pi/2,np.pi/2,len(strain))
    win=np.cos(x)**1


    # FT of data
    noise_ft=np.fft.fft(win*strain)

    # Smoothing our data to create weights
    noise_smooth=smooth_vector(np.abs(noise_ft)**2, 10)
    noise_smooth=noise_smooth[:len(noise_ft)//2+1] #will give us same length
    tobs=dt*len(strain)
    dnu=1/tobs
    nu=np.arange(len(noise_smooth))*dnu
    nu[0]=0.5*nu[1]
    Ninv=1/noise_smooth     # Creating weights
    Ninv[nu>1500]=0
    Ninv[nu<20]=0

    # Whitening the data
    whitening_range = np.arange(1100, 5000)
    vs_noise = smooth_vector(np.abs(noise_ft)**2, 10)[whitening_range]
    w_noise = noise_ft[whitening_range]/vs_noise
    plt.loglog(np.abs(noise_ft)**2)
    plt.loglog(vs_noise)
    plt.loglog(np.abs(w_noise)**2)
    
    

    #plt.show()

    # FFT of template
    template_ft=np.fft.rfft(tp*win)
    template_filt=template_ft*Ninv

    # Matched filtering
    data_ft=np.fft.rfft(strain*win)

    filter = np.ones(len(data_ft))
    filter[nu > 500] = 0
    filter[nu < 40] = 0
    #filter = filter*0 + 1
    plt.clf()
    plt.loglog(nu, np.abs(data_ft*filter)**2)
    plt.show()


    #rhs=np.fft.irfft(data_ft*np.conj(template_filt)*filter)
    rhs=np.fft.irfft(data_ft*np.conj(template_filt))

    return rhs

#h = getnobel('LIGO/H-H1_LOSC_4_V2-1135136334-32.hdf5')
l = getnobel('LIGO/L-L1_LOSC_4_V2-1135136334-32.hdf5', tmp_name = "LIGO/GW170104_4_template.hdf5")

plt.clf()
#plt.plot(h)
plt.plot(l, ":k")
plt.show()

