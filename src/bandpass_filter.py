import wave, struct
import numpy as np
from scipy.signal import butter, lfilter
import matplotlib.pyplot as plt
import sys
from optparse import OptionParser
import os

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

def writefile(fname, fileparams, data):
    wav_file = wave.open(fname, "w")
    wav_file.setparams(fileparams)
    for d in data:
        wav_file.writeframes(struct.pack('h',d))
    wav_file.close()



def main(argv):

   # input parser
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.set_defaults(order=6)
    parser.set_defaults(low_cut=1e3)
    parser.set_defaults(high_cut=50e3)
    parser.set_defaults(Fs=144000)

    parser.add_option("-s", "--sampling", type="float", dest="Fs",
                      help="sampling frequency, default = 144000")
    parser.add_option("-o", "--order", type="int", dest="order",
                      help="Order of Butterwordth, default = 6")
    parser.add_option("-f", "--low", type="float", dest="low_cut",
                      help="Lower cut off  [Hz], default = 1e3")
    parser.add_option("-F", "--high", type="float", dest="high_cut",
                      help="High cut off, default = 50e3")

    (options, args) = parser.parse_args()


    # Filter requirements.
    Fs = options.Fs
    order = options.order
    cutoff_low = options.low_cut  # desired cutoff frequency of the filter, Hz
    cutoff_high = options.high_cut   # desired cutoff frequency of the filter, Hz
    
    waveFile = wave.open('myfile2.wav', 'r')
    length = waveFile.getnframes()
    time_series = []

    for j in range(0, length/10000):
        waveData = waveFile.readframes(1)
        sample_point = struct.unpack("<h", waveData)
        time_series.append(sample_point[0])


    # Filter the data, and plot both the original and filtered signals.
    time_series_filtered = butter_bandpass_filter(time_series, cutoff_low, cutoff_high, Fs, order)

    path = sys.argv[-1]
    basename = os.path.splitext(os.path.basename(path))[0]
    writefile(basename+'_filtered.wav', waveFile.getparams(), time_series)

    
if __name__ == "__main__":
    main(sys.argv)


