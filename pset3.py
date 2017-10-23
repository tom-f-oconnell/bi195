#!/usr/bin/env python

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('dark')

df = pd.read_csv('data_pset3.csv')

data = df['x']
times = df['secs']
dt = times[1] - times[0]

power = np.abs(np.fft.fft(data))**2

# TODO what exactly is fftfreq doing?

# d=sample spacing, specific units irrelevant
freqs = np.fft.fftfreq(len(data), d=dt)
nonnegative_indices = freqs > 0

plt.semilogy(freqs[nonnegative_indices], power[nonnegative_indices])
plt.title('Sound power in animal facility')
plt.xlabel('Frequency ($Hz$)')
plt.ylabel('Power')
plt.savefig('sound_power_pset3.svg')
plt.show()
