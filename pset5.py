#!/usr/bin/env python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('dark')

n_repeats = 1000
n_draws = 100

xs = np.zeros(n_repeats) * np.nan

for i in range(n_repeats):
    # [0,1)
    draws = np.random.rand(n_draws)
    xs[i] = np.sum(draws)

assert np.sum(np.isnan(xs)) == 0

# the mean should be n_draws * E[single_draw] = 50
# the standard deviation should be 
