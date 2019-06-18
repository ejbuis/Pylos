from math import *
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from scipy import signal
import sys


def cosmic_ray(amplitude=0.05):
    '''
    EJB: generate sound that is induced by a cosmic ray particle in water
    '''
    freq = 20000.
    pulseheight = 0.
    pulsewidth = 10
    for i in count(0):
        x =  (i%freq)-int(freq/2)
        iks = float(x)
        norm = pulsewidth*math.exp(-1.)
        if x > 10000:
            pulseheight = 0.0
        else:
            aid = -(iks*iks)/(pulsewidth*pulsewidth)
            pulseheight = -iks*math.exp(aid)
#        print x, ((pulseheight)/norm) * 0.01#amplitude
        yield ((pulseheight)/norm) * 0.75#amplitude


def main(argv):
    
    sampling_rate = 100000
    time_series = clicks[1]
    t = np.arange(len(time_series))
    f, spectrum = PSD_using_scipy(time_series, sampling_rate)

    fig, (ax0, ax1) = plt.subplots(2, 1)
    ax0.plot(t, time_series)
    ax1.plot(f, spectrum) 
    ax0.grid(True)
    ax1.grid(True)
    ax0.set_xlabel('sample nr')
    ax1.set_xlabel('frequency [Hz]')
    ax0.set_ylabel('ADC')
    ax1.set_ylabel('PSD [ADC/$\sqrt{Hz}$]')
    plt.subplots_adjust(hspace=0.4)
    plt.show()

        from scipy.io.wavfile import write
    clicks = write("neutrinos.wav")



if __name__ == '__main__':
    main(sys.argv)

    
