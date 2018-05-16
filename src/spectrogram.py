import wave, struct
from math import *
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from scipy import signal
import sys
import pylab
import os


def graph_spectrogram(wav_file, begin, end):
    sound_info, Fs = get_wav_info(wav_file)
    pylab.figure(num=None, figsize=(19, 12))
    pylab.subplot(111)
    pylab.title('spectrogram of %r' % wav_file)
    print len (sound_info[begin*Fs: end*Fs]), Fs, len(sound_info)
    pylab.specgram(sound_info[begin*Fs: end*Fs], Fs=Fs)
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
    end = 10000
    if len(argv) > 1:
        filename  = argv[1]
    if len(argv) > 2:
        begin = string.atoi(argv[2]) # in seconds
        end = string.atoi(argv[3]) # seconds
    graph_spectrogram(argv[1], begin, end)

if __name__ == '__main__':
    main(sys.argv)
