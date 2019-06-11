import wave, struct
from math import *
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from scipy import signal
import sys
import pylab
import os

import datetime
import string

from bandpass_filter import *

def usage():
    print("Usage:  %s  <wav-file> <t_min (s)> <t_max (s)> \n"
          % os.path.basename(sys.argv[0]))


def get_timestamp(filename, begin_time):
    s = os.path.splitext(os.path.basename(filename))[0][10:]
    extra_minute = int(begin_time/60.)
    extra_second = int(begin_time - extra_minute*60.)
    hours = int(s[6:8])
    minutes = int(s[8:10]) + extra_minute
    seconds = int(s[10:12]) + extra_second
    if seconds > 60:
        minutes += 1
        seconds -= 60
    if minutes > 60:
        hours +=1
        minutes -= 60
        
    date = datetime.datetime(year=2000+ int(s[0:2]), month=int(s[2:4]), \
                             day=int(s[4:6]), \
                             hour = hours, minute = minutes, second = seconds)
    return date

def graph_spectrogram(wav_file, begin, end, cutoff_low, cutoff_high, order):
    sound_info, Fs = get_wav_info(wav_file)
    pylab.figure(num=None, figsize=(19, 12))
    pylab.subplot(111)
    pylab.title('spectrogram of %r' % wav_file)
    nfft = 1024
    data = sound_info[int(begin*Fs): int(end*Fs)]
    data = butter_bandpass_filter(data, cutoff_low, cutoff_high, Fs, order)
    pylab.specgram(data, Fs=Fs, NFFT = nfft, \
                   noverlap =int(nfft/2),
                   cmap='nipy_spectral', \
                   mode = 'psd')
    pylab.colorbar()
    
    #    datestring = get_timestamp(wav_file, begin)
    #    pylab.xlabel("Time (seconds since {}) [s]".format(datestring), \
        #                 size = 16, ha ='right', x=1.0)
    pylab.xlabel("Time (seconds)", size = 16, ha ='right', x=1.0)
    pylab.ylabel("PSD", size = 16, ha ='right', position=(0,1))
    pylab.savefig('spectrogram.png')
    pylab.show()
    
def get_wav_info(wav_file):
    wav = wave.open(wav_file, 'r')
    frames = wav.readframes(-1)
    sound_info = pylab.fromstring(frames, 'int16')
    frame_rate = wav.getframerate()
    wav.close()
    return sound_info, frame_rate


import string
def main(argv):
    begin = 0
    end = 100
    
    if len(argv) < 2:
        usage()
        sys.exit()
    if len(argv) > 1:
        filename  = argv[1]
    if len(argv) > 2:
        begin = float(argv[2]) # in seconds
        end = float(argv[3]) # seconds

    filename = argv[1]
    cutoff_low = 100.
    cutoff_high = 30000.
    order = 5 
    graph_spectrogram(filename, begin, end, cutoff_low, cutoff_high, order)

if __name__ == '__main__':
    main(sys.argv)
