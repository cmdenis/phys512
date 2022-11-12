import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
import json


def smooth_vector(vec,sig):
    n = len(vec)                               
    x = np.arange(n)                          
    x[n//2:] = x[n//2:]-n                       # Axis
    kernel = np.exp(-0.5*x**2/sig**2)           # Gaussian to convolve
    kernel = kernel/kernel.sum()                # Normalize Gaussian
    vecft = np.fft.rfft(vec)                    # FFT of array
    kernelft = np.fft.rfft(kernel)              # FFT of gaussian
    vec_smooth = np.fft.irfft(vecft*kernelft)   # convolve the data with the kernel
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



def getlocation(ename):

    print("Looking at event", ename)

    # Opening JSON file to map and getting info
    with open("LIGO/BBH_events_v3.json") as f:
        events = f.read()

    events = json.loads(events)
    event = events[ename]

        
    tmp_name = "LIGO/" + event["fn_template"]
    fname_h = "LIGO/" + event["fn_H1"]
    fname_l = "LIGO/" + event["fn_L1"]
    fname_list = [fname_h, fname_l]

    locations = []
    snr = []
    ana_snrs = []

    for fname in fname_list:
        # Loading the data
        print('Reading File:',fname)
        strain,dt,utc=read_file(fname)

        # Loading the template
        print("Using template:", tmp_name)
        tp,tx=read_template(tmp_name)
        
        # Generating the cosine window
        x=np.linspace(-np.pi/2,np.pi/2,len(strain))
        win=np.cos(x)

        # Taking the fft of the data
        noise_ft=np.fft.fft(win*strain)

        

        # Creating a smooth version of the spectrum
        noise_smooth=smooth_vector(np.abs(noise_ft)**2, 20)

        # Whitened data
        w_data = noise_ft/noise_smooth

        # Generating the axis
        tobs=dt*len(strain)
        dnu=1/tobs
        nu=np.arange(len(noise_smooth))*dnu
        nu[0]=0.5*nu[1]

        w_data[nu < 10] = 0
        w_data[nu > 1000] = 0
        

        # Match filtering
        template_ft=np.fft.fft(tp*win)
        rhs=np.fft.irfft(w_data*np.conj(template_ft))
        times = np.arange(len(rhs))*dt/2


 
        # Analytical SNR
        ana_ns = 1/np.sqrt(np.sum(np.abs(template_ft)**2/noise_smooth))
        ana_snr = np.mean(np.abs(template_ft)**2)/ana_ns


        # Calculating the SNR
        window = 2500
        int_range = rhs[np.argmax(np.abs(rhs))-window:np.argmax(np.abs(rhs))+window]
        int_times = times[np.argmax(np.abs(rhs))-window:np.argmax(np.abs(rhs))+window]
        



        noise_amp = np.mean(int_range**2)
        signal_amp = int_times[np.argmax(np.abs(int_range))]**2

        snr.append(signal_amp/noise_amp)
        ana_snrs.append(ana_snr)

        locations.append(times[np.argmax(np.abs(rhs))])

    return locations, snr, ana_snrs

event_list = ['GW150914', "LVT151012", "GW151226", "GW170104"]

for event in event_list:    
    print()
    times, snr, ana_snr = getlocation(event)
    print("Time of events at both detectors:", times, "with difference: ", np.abs(times[1] - times[0]))
    print("SNRs:", snr)
    print("Combined SNR:", ((np.sqrt(snr[0])+np.sqrt(snr[1]))**2)/2)
    print("Analytical SNR:", ana_snr)
