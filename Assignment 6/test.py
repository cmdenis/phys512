import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
import json


def smooth_vector(vec,sig):
    n=len(vec)
    x=np.arange(n)
    x[n//2:]=x[n//2:]-n
    kernel=np.exp(-0.5*x**2/sig**2) # make a Gaussian kernel
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
        w_data[nu > 800] = 0
        

        # Plotting
        plt.loglog(nu, np.abs(noise_ft)**2)
        plt.loglog(nu, noise_smooth)
        plt.loglog(nu, np.abs(w_data)**2)
        plt.show()

        # Match filtering
        template_ft=np.fft.fft(tp*win)
        rhs=np.fft.irfft(w_data*np.conj(template_ft))
        times = np.arange(len(rhs))*dt/2

        # Plotting
        plt.plot(times, rhs)
        plt.xlabel("Time (s)")
        plt.ylabel("m")
        plt.savefig("figs/a6q5_example_location.jpg")
        plt.show()

        locations.append(times[np.argmax(np.abs(rhs))])


    return locations

event_list = ['GW150914', "LVT151012", "GW151226", "GW170104"]

for event in event_list:    
    times = getlocation(event)
    print("Time of events at both detectors:", times)
    print("Difference in time", np.abs(times[1] - times[0]))


assert(0==1)


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