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
    gpsStart=meta['GPSstart'][()]
    utc=meta['UTCstart'][()]
    duration=meta['Duration'][()]
    strain=dataFile['strain']['Strain'][()]
    dt=(1.0*duration)/len(strain)
    dataFile.close()
    return strain,dt,utc



def locate_peak(ename, tmp_name = 'LIGO/GW150914_4_template.hdf5'):

    # Opening JSON file to map and getting info
    with open("LIGO/BBH_events_v3.json") as f:
        events = f.read()

    events = json.loads(events)
    event = events[ename]

    
    tmp_name = "LIGO/" + event["fn_template"]
    fname_h = "LIGO/" + event["fn_H1"]
    fname_l = "LIGO/" + event["fn_L1"]
    fname_list = [fname_h, fname_l]

    for fname in fname_list:
        # Loading data
        print('reading file ', fname)
        strain, dt, utc = read_file(fname)

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
        vs_noise = smooth_vector(np.abs(noise_ft)**2, 300)[whitening_range]
        w_noise = noise_ft[whitening_range]/vs_noise
        


        # FFT of template
        template_ft=np.fft.rfft(tp*win)
        template_filt=template_ft*Ninv

        # Matched filtering
        data_ft=np.fft.rfft(strain*win)

        filter = np.arange(len(data_ft))
        filter[filter < 1000] = 0
        filter[filter > 10000] = 0
        filter[filter > 0] = 1
        filter = filter*0 + 1
    


        #rhs=np.fft.irfft(data_ft*np.conj(template_filt)*filter)
        rhs=np.fft.irfft(data_ft*np.conj(template_filt))

        peak_location = np.argmax(np.abs(rhs))

    return peak_location, dt, len(rhs)


# Pair events list:
events = [
    ['LIGO/H-H1_LOSC_4_V1-1167559920-32.hdf5', 'LIGO/L-L1_LOSC_4_V1-1167559920-32.hdf5'],
    ['LIGO/H-H1_LOSC_4_V2-1126259446-32.hdf5', 'LIGO/L-L1_LOSC_4_V2-1126259446-32.hdf5'],
    ['LIGO/H-H1_LOSC_4_V2-1128678884-32.hdf5', 'LIGO/L-L1_LOSC_4_V2-1128678884-32.hdf5'],
    ['LIGO/H-H1_LOSC_4_V2-1135136334-32.hdf5', 'LIGO/L-L1_LOSC_4_V2-1135136334-32.hdf5']
]


for i in events:
    h, dt, length = locate_peak(i[0])
    l, dt, length = locate_peak(i[1])

    diff1 = np.abs(h-l)
    diff2 = np.abs(np.abs(d))

    print("Hanford:", h)
    print("Livingston:", l)
    print("Time difference", (h-l)*dt)

