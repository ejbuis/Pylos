import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as fft
from scipy import signal

def PSD_using_numpy(ampl, fs):
    print 'sample freq ', fs
    spectrum = fft.fft(ampl)
    PSD = np.abs(spectrum)**2 * (2./(fs*len(spectrum)))
    freq = fft.fftfreq(len(spectrum),d=1./fs)
    return freq, PSD

def main():
    bip = np.array([ -1.7513e-37, -3.37337e-35, -4.04037e-41, 2.50476e-30, 9.65528e-28, 2.21075e-25, 3.54893e-23, 4.18286e-21, 3.66047e-19, 2.3519e-17, 1.05946e-15, 2.83482e-14, 3.14432e-21, -4.06164e-11, -2.17489e-09, -6.91748e-08, -1.54256e-06, -2.52555e-05, -0.000307013, -0.00274015, -0.0171466, -0.0637317, -4.90981e-10, 1.762, 13.1062, 57.9064, 179.374, 407.952, 688.885, 854.087, 742.406, 383.316, -0, -204.494, -211.296, -129.681, -55.8014, -17.6292, -4.1353, -0.712197, -0.0859956, -0.00616778, 9.1688e-13, 6.34932e-05, 9.11326e-06, 7.76956e-07, 4.64411e-08, 2.03811e-09, 6.64108e-11, 1.5888e-12, 2.6649e-14, 2.65504e-16, -1.09653e-26, -5.27405e-20, -1.05154e-21, -1.24534e-23, -1.03402e-25, -6.30365e-28, -2.85325e-30, -9.48216e-33, -2.20931e-35, -3.05762e-38, 2.63126e-49, 1.17201e-43 ])

    Fs = 1/(4e-6)
    t = np.linspace(0, len(bip)/Fs, len(bip))
    freq, freq_data = PSD_using_numpy(bip, Fs)

    Fs_resampled = 144000.
    number_resampled = int(round(len(bip)*Fs_resampled/Fs))
    bip_resampled = signal.resample(bip, number_resampled)
    t_resampled = np.linspace(0, len(bip_resampled)/Fs_resampled, len(bip_resampled)) 
    freq_resampled, freq_data_resampled = PSD_using_numpy(bip_resampled, Fs_resampled)
    
    plt.subplot(411)
#    plt.plot(t, bip)
    plt.plot(t_resampled, bip_resampled, '--', c='blue')
    plt.plot(bip_resampled, '--', c='blue')

    plt.subplot(412)
    plt.plot(freq, freq_data)
    plt.plot(freq_resampled, freq_data_resampled, '--', c='blue')

    sig = np.tile(bip_resampled, 10)
    sig_noise = sig + 0.1* bip_resampled.max()* np.random.randn(len(sig)) #dit moet echte data zijn...
    plt.subplot(413)
    plt.plot(sig_noise)

    fir_coefficients = bip_resampled[22:12:-1]
    det = signal.lfilter(fir_coefficients, 1, sig_noise)
    det = det*det
    plt.subplot(414)
    plt.plot(det)
    plt.show()

import sys
if __name__ == "__main__":
    sys.exit(main())

