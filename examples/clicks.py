from math import *
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from scipy import signal
import sys

def PSD_using_scipy(ampl, fs):
    return signal.welch(ampl, fs, nperseg=1024, scaling= 'spectrum')

def main(argv):
    from scipy.io.wavfile import read
    clicks = read("spermwhale.wav")
    
    sampling_rate = clicks[0]
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


if __name__ == '__main__':
    main(sys.argv)
