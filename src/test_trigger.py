import wave, struct
from math import *
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from scipy import signal
import sys

def writefile(fname, fileparams, data):
    wav_file = wave.open(fname, "w")
    wav_file.setparams(fileparams)
    for d in data:
        wav_file.writeframes(struct.pack('h',d))
    wav_file.close()

def PSD_using_scipy(ampl, fs):
    return signal.welch(ampl, fs, nperseg=1024, scaling= 'spectrum')

def main(argv):
    waveFile = wave.open('myfile2.wav', 'r')
    length = waveFile.getnframes()
    time_series = []
    print   waveFile.getparams()

    for i in range(0,95):
        waveFile.setpos(int(i*length/100))
        for j in range(0, (length/100)):
            waveData = waveFile.readframes(1)
            sample_point = struct.unpack("<h", waveData)
            if sample_point[0] > 15000:
                this_point = waveFile.tell()
                print "Trigger ", i, this_point, sample_point[0]
                waveFile.setpos(this_point - 20000)
                for k in range(0, 30000):
                    waveData = waveFile.readframes(1)
                    sample_point = struct.unpack("<h", waveData)
                    time_series.append(sample_point[0])

    writefile("test.wav", waveFile.getparams(), time_series)
    sampling_rate = waveFile.getframerate()
    t = np.arange(len(time_series))

    f, spectrum = PSD_using_scipy(time_series, sampling_rate)

    fig, (ax0, ax1) = plt.subplots(2, 1)
    ax0.plot(t, time_series)
    ax1.loglog(f, spectrum) 
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
