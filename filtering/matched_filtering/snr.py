import matplotlib.pyplot as pp
    
import pycbc.filter
import pycbc.psd
from pycbc import types
from pycbc.types.array import complex_same_precision_as
from pycbc import fft
import sys

import numpy as np
from numpy.linalg import norm as Norm
import wave

from scipy import signal
from scipy.io.wavfile import read
from scipy.signal import butter, lfilter

import math

Nperseg = int(2048)

def butter_highpass(cutoff, Fs, order=5):
    nyq = 0.5 * Fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, Fs, order=5):
    b, a = butter_highpass(cutoff, Fs, order=order)
    y = lfilter(b, a, data)
    return y

def Knudsen_curve(seastate = 0, frequency=1000):
    assert seastate < 7
    ss = [44.5, 50, 55., 61.5, 64.5, 66.5, 70]
    i = 1
    if frequency >0 :
        return ss[seastate] -17.* math.log10(frequency/1000.)
    else:
        return 1

def frequency_span(Fs, NFFT):
    return np.linspace(0, int(NFFT/2-1), int(NFFT/2))*Fs/NFFT

def spectrum(NFFT, Fs, data):
    length     = len(data)
    PSD        = np.zeros((int(NFFT/2), 1), dtype=float)
    freq       = frequency_span(Fs, NFFT)
    Segment    = int(length/NFFT)
    print("segment length =", NFFT)
    wnd        = np.hanning(NFFT)
    norm       = Norm(wnd)**2
    double2single   = 2.0
    for span in range(0, Segment):
        Bg          = span*NFFT
        end         = Bg+NFFT
        yw          = wnd*data[Bg:end]
        a           = np.fft.fft(yw, NFFT)
        ac          = np.conj(a)
        pxx         = np.abs(a*ac)
        PSD[:, 0]  +=  double2single*pxx[0:int(NFFT/2)]
        PSD[:, 0]  /= (float(Segment)*NFFT*norm) 
        #return 10*np.log10(PSD[:, 0])
        return types.FrequencySeries(10*np.log10(PSD[:, 0]), freq[1] - freq[0])

    
def datafile(filename):
    # read one column time trace (sample freq 144 kHz)
    data = np.loadtxt(filename,  usecols=(0), dtype='float', unpack=True)  
    dt = 1/144000.
    data_td = types.TimeSeries(data, dt)

    return data_td

def padme(array1, array2):
    # make length as the longest array
    if len(array1) > len(array2):
        array2 = np.pad(array2, (len(array1)-len(array2), 0) )
        return array1, array2
    elif len(array2) > len(array1):
        array1 = np.pad(array1, (len(array2)-len(array1), 0) )
        return array1, array2
    else:
        return array1, array2

def main(argv):
    # quick and dirty script to get going...
    # Need to get a SNR array using matched filtering
    # Matched filtering: template, data and a noise PSD

    template_file = 'neutrino_template.dat'
    data_file = 'neutrino_6_300_71.txt'    
    template_array = datafile(template_file)
    data_array = datafile(data_file) 

    # take care of file length
    #    template_array, data_array = padme(template_array, data_array)

    # convert
    fs = 144000.
    dt = 1/fs
    data = types.TimeSeries(data_array, dt)
    template = types.TimeSeries(template_array, dt)
    template.resize(len(data))

    # high pass filter for the data
    filter_cut = 100
#    time_series_filtered = butter_highpass_filter(data, filter_cut, fs ,order = 5) 

    ### scalings factors
    scaling_factor_min = 1e3
    scaling_factor_max = 2e4
    psd3 = spectrum(Nperseg, fs, data*scaling_factor_min)
    psd  = spectrum(Nperseg, fs, data*scaling_factor_max)
#
#    # plot time trace
#    ax2 = pp.gca()
#    
#    ax2.plot(data.sample_times, time_series_filtered*scaling_factor_min, label='data time trace')
#    ax2.set_ylabel('pressure ($\mu$Pa)')
#    ax2.set_xlabel('Time (s)')
#    pp.show()
#    pp.close()
#
#
#    # Knudsen curves: unit is db Re 1muPa^2/Hz
#    Knudsen_noise_0 = np.array([Knudsen_curve(0,i) for i in psd.sample_frequencies])
#    Knudsen_noise_6 = np.array([Knudsen_curve(6,i) for i in psd.sample_frequencies])
#
#    #
    ax = pp.gca()
#    
#    ax.semilogx(psd.sample_frequencies, psd, label='psd')
#    ax.semilogx(psd3.sample_frequencies, psd3, label='psd')
#
#    ax.semilogx(psd.sample_frequencies, Knudsen_noise_0, label='sea state 0')
#    ax.semilogx(psd.sample_frequencies, Knudsen_noise_6, label='sea state 6')
#
#    pp.xlim(1000,100000)
#    ax.legend()
#    ax.grid()
#    ax.set_ylabel('PSD (dB Re 1$\mu$Pa$^2$/Hz)')
#    ax.set_xlabel('Frequency (Hz)')
#

    p_estimated = data.psd( Nperseg*1/fs, avg_method='mean',  window='hann') 
    
    p = pycbc.psd.interpolate(p_estimated, data.delta_f)
#    p = psd.inverse_spectrum_truncation(p, int(dt*seg_len*noise_td.sample_rate),low_frequency_cutoff=flow)
#    psd_inter = p

    if 'match' in argv:
        flow = 1000
        fhigh = 25000
        snr = pycbc.filter.matched_filter(template*i, data,
                                        psd = p,
                                          low_frequency_cutoff=flow, high_frequency_cutoff=fhigh)

        ax.plot(snr.sample_times, abs(snr), label='snr')
    pp.show()
    pp.close()
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))


