#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 14:07:51 2018

@author: rosavandenende
"""

#in dit bestand: sigma bepalen, en vervolgens kijken welke ADC waarden boven 5*sigma liggen. 


import wave, struct
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import butter, lfilter
import glob 

def butter_lowpass(cutoff, sampling_rate, order=5):
    nyq = 0.5 * sampling_rate
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, sampling_rate, order=5):
    b, a = butter_lowpass(cutoff, sampling_rate, order=order)
    y = lfilter(b, a, data)
    return y

def butter_highpass(cutoff, sampling_rate, order=5):
    nyq = 0.5 * sampling_rate
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, sampling_rate, order=5):
    b, a = butter_highpass(cutoff, sampling_rate, order=order)
    y = lfilter(b, a, data)
    return y

def butter_bandpass(lowcut, highcut, sampling_rate, order=5):
    nyq = 0.5 * sampling_rate
    low = lowcut/nyq
    high = highcut/nyq
    print (low, high)
    b, a = butter(order, [low, high], btype='band', analog=False)
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, sampling_rate, order=5):
    y = butter_lowpass_filter(data, highcut, sampling_rate, order)
    return butter_highpass_filter(y, lowcut, sampling_rate, order)

def PSD_using_scipy(ampl, fs):
    return signal.welch(ampl, fs, nperseg=1024, scaling= 'spectrum')

def main_trigger(n,a,b,low_cut,high_cut):
    file_list = glob.glob("/Volumes/My Passport/Pylos_data/Py1-201650198/*.wav")                    #Maakt een lijst met alle bestandsnamen uit het mapje
    time_series = []
    samples_per_sigma = 120 * 144000                                                                #Elke 120 seconden wordt sigma opnieuw bepaald
    sigma_list = []
    sigma2_list = []
    triggered = 0
    triggered_list = []
    teller = 0 
    teller_list = []
   
    for i in file_list[a:b]:
        waveFile = wave.open(i,'r')
        sampling_rate = waveFile.getframerate()
        for j in range(0,int(n)):                                  
            waveFile.setpos(j)
            waveData = waveFile.readframes(1)
            sample_point = struct.unpack("<h", waveData)
            time_series.append(sample_point[0]) 

    time_series_filtered = butter_bandpass_filter(time_series,low_cut,high_cut, sampling_rate,order = 5) 

    for i in range(0,int(n)//int(samples_per_sigma)):                                               #Grootste deler                              
        sigma = np.std(time_series_filtered[samples_per_sigma*i : samples_per_sigma*(i+1)])         #Bepaling van sigma voor elk stuk van 120 sec       
        sigma_list.append(sigma)
        
        for element in time_series_filtered[samples_per_sigma*i : samples_per_sigma*(i+1)]:         #Voor elk element in dit stuk lijst:  trigger?
            teller += 1                                                                             #Teller: op welk sample number bevindt de trigger zich 
            if element >= (5*sigma):
                triggered += 1
                triggered_list.append(element)
                teller_list.append(teller) 
    
    #als niet deelbaar door samples_per_sigma, dan moet je nog het stuk hebben van k+1 * samples_per_sigma tot het einde, waarbij k de grootste deler is 
    rest = n - (int(n)//int(samples_per_sigma)+1)*samples_per_sigma                                                        
    if rest != 0:
        sigma = np.std(time_series_filtered[rest:])
        sigma_list.append(sigma)
        rest_len = len(time_series_filtered[rest:])
        for element in time_series_filtered[rest:]:
            teller += 1
            if element >= (5*sigma):
                triggered += 1
                triggered_list.append(element)
                teller_list.append(teller) 

    for sigma_value in sigma_list[:-1]:                                                               #Een lijst maken met de sigmavalues, evenvaak als de samples
        sigma2_list.extend([5*sigma_value]*(samples_per_sigma)) 
    sigma2_list.extend([5*sigma_list[-1]]*(rest_len))
     
    plt.rcParams["figure.figsize"] = [24,12]
    fig,ax = plt.subplots(2,1,sharex = False)
    ax[0].scatter(teller_list,triggered_list,s=10)
    ax[0].set_xlabel("sample number")
    ax[0].set_ylabel("ADC")
    ax[0].set_title("Triggered signals, amplitude > 5 * $\sigma$")
    ax[0].text(100000,350, "Number of triggers = %d" %triggered)
    
    ax[1].plot(range(len(time_series_filtered)),time_series_filtered,'green',linewidth = 0.7)
    ax[1].plot(range(len(time_series_filtered)),sigma2_list,'r--')
    ax[1].set_title("Filtered signal (1-40 kHz)")
    ax[1].set_xlabel("sample number")
    ax[1].set_ylabel("ADC")
    plt.show()

main_trigger(3*60*144000,3,4,1e3,40e3)
