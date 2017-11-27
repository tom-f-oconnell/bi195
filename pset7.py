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

def autocorr_and_welch(data, segment_length=None):
    if segment_length is None:
        f1 = plt.figure()
        # "compute the correlation function of this signal and plot it"
        autocorr = np.correlate(data, data, mode='same')

        plt.title('Autocorrelation as a function of time lag')
        plt.xlabel('Time lag')
        plt.ylabel('Correlation')
        plt.plot(autocorr)

    # welch
    if not segment_length is None:
        freqs, powers = scipy.signal.welch(data, nperseg=segment_length)
    else:
        freqs, powers = scipy.signal.welch(data)

    f2 = plt.figure()
    plt.loglog(freqs, powers)
    plt.xlabel('Frequency')
    plt.ylabel('Power')
    title = 'PSD estimate'
    if not segment_length is None:
        title += ' segment_length={}'.format(segment_length)
    plt.title(title)

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

autocorr_and_welch(summed, segment_length=1024)
autocorr_and_welch(summed, segment_length=64)

plt.show()
