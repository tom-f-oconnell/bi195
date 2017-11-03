#!/usr/bin/env python

from __future__ import division
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('dark')

"""
Problem 3
"""
n_repeats = 1000
n_draws = 100

xs = np.zeros(n_repeats) * np.nan

for i in range(n_repeats):
    # [0,1)
    draws = np.random.rand(n_draws)
    xs[i] = np.sum(draws)

assert np.sum(np.isnan(xs)) == 0

# the mean should be n_draws * E[single_draw] = 50
# the standard deviation should be sqrt(n_draws / 12)

manual_mean = n_draws / 2
# 5 / np.sqrt(3) in this case
manual_stddev = np.sqrt(n_draws / 12.0)

numeric_mean = np.mean(xs)
numeric_stddev = np.std(xs)

print 'Numeric mean: {}'.format(numeric_mean)
print 'Numeric stddev: {}'.format(numeric_stddev)
print 'Manual mean: {}'.format(manual_mean)
print 'Manual stddev: {}'.format(manual_stddev)
print 5 / np.sqrt(3)

# mean
mu = manual_mean
sigma = manual_stddev

sigmas_to_show = 4
xx = np.linspace(mu - sigmas_to_show * sigma, mu + sigmas_to_show * sigma, num=1000)
plt.plot(xx, norm.pdf(xx, mu, sigma), 'r', linewidth=2, alpha=0.9)

bins = 50
n, _, _ = plt.hist(xs, bins, normed=1, facecolor='black', alpha=0.75)
plt.ylim([0, max(n) * 1.2])
plt.title('CLT: Normal vs. empirical distribution')
plt.xlabel('$x$')
plt.ylabel('Frequency')

plt.show()
