import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as fft
from scipy import signal
from scipy.signal import butter, lfilter

def PSD_using_scipy(ampl, fs):
    return signal.welch(ampl, fs, nperseg=1024, scaling= 'spectrum')

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut/nyq
    high = highcut/nyq
    print low, high
    b, a = butter(order, [low, high], btype='band', analog=False)
    return b, a

#def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
#    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
#    y = lfilter(b, a, data)
#    return y

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    y = butter_lowpass_filter(data, highcut, fs, order)
    return butter_highpass_filter(y, lowcut, fs, order)

def PSD_using_numpy(ampl, fs):
    spectrum = fft.fft(ampl)
    PSD = np.abs(spectrum)**2 * (2./(fs*len(spectrum)))
    freq = fft.fftfreq(len(spectrum),d=1./fs)
    return freq, PSD

def main():
    Fs = 144000
    
    template = np.loadtxt('template_click.dat')
    data = np.loadtxt('sample-data2.txt')
    sig_noise = butter_highpass_filter(data, 1000 , Fs, 5)
    
    fir_coefficients = template[1250:1200:-1]
    det = signal.lfilter(fir_coefficients, 1, sig_noise)

    f, det_spec = PSD_using_scipy(det*det, 1)
    
    plt.subplot(311)
    plt.plot(sig_noise)
    plt.subplot(312)
    plt.plot(det*det)
    plt.subplot(313)
    plt.plot(f, det_spec, 'o-')
    plt.xlim([0,0.1])

#    plt.hist(det, 50)
    plt.show()

import sys
if __name__ == "__main__":
    sys.exit(main())

