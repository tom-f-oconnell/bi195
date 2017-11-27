#!/usr/bin/env python

# Bi195 - Mathematics in Biology
# Problem set 7
# Thomas O'Connell

from __future__ import division
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('dark')

# Question 1: Power spectrum of a stationary Gaussian process

def autocorrelation(data):
    num_valid_lags = len(data) - 1
    autocorr = np.zeros(num_valid_lags)
    for i in range(1, num_valid_lags):
        '''
        print data[:-i].shape
        print data[i:].shape
        '''
        autocorr[i] = np.corrcoef(data[:-i], data[i:])[0, 1]
    return autocorr

def autocorr_and_welch(data):
    plt.figure()
    plt.plot(data)

    f1 = plt.figure()
    # "compute the correlation function of this signal and plot it"
    #autocorr = np.correlate(data, data, mode='same')

    autocorr = autocorrelation(data)

    plt.title('Autocorrelation as a function of time lag $\tau$)')
    plt.xlabel('$\tau$')
    plt.ylabel('Correlation')
    plt.plot(autocorr)

    # welch
    f2 = plt.figure()
    freqs, powers = scipy.signal.welch(data)
    plt.loglog(freqs, powers)
    plt.xlabel('Frequency')
    plt.ylabel('Power')
    plt.title('PSD estimate')

samples = 100000
# from mean=0, variance=std=1
xt = np.random.randn(samples)

autocorr_and_welch(xt)

bin_size = 100
summed = np.zeros(samples - bin_size + 1) * np.nan
for t in range(len(summed)):
    assert len(xt[t:t+100]) == 100
    summed[t] = np.sum(xt[t:t+100])

assert np.sum(np.isnan(summed)) == 0

autocorr_and_welch(summed)

# TODO do welch bonus as well

#autocorr = np.correlate(xt, xt, mode='full')
#f2 = plt.figure()
#plt.plot(autocorr)

plt.show()
