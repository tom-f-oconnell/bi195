#!/usr/bin/env python

from __future__ import division
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('dark')

"""
Problem 3: CLT
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
'''
plt.plot(xx, norm.pdf(xx, mu, sigma), 'r', linewidth=2, alpha=0.9)

bins = 50
n, _, _ = plt.hist(xs, bins, normed=1, facecolor='black', alpha=0.75)
plt.ylim([0, max(n) * 1.2])
plt.title('CLT: Normal vs. empirical distribution')
plt.xlabel('$x$')
plt.ylabel('Frequency')

plt.show()
'''

"""
Problem 4: PCA
"""

# samples / measurements
m = 1000

# covariance diagonal
# since these covariance matrices are zero off the diagonal
# the r.v.s of each component are independent
# TODO is that true? we could still have marginal normals and zero covariances,
# but a non-normal joint distribution where the r.v.s are dependent, right?
C_ii = [1, 0.25]
n = len(C_ii)

# TODO is this the usual convention for dimensions? (which is measurement #)
M = np.zeros((m,n)) * np.nan

for i in range(n):
    # to get samples from N(mu,sigma^2) via the "standard normal" (i.e. N(0,1))
    # sigma * np.random.randn(...) + mu
    M[:,i] = np.sqrt(C_ii[i]) * np.random.randn(m)

# computing the empirical covariance matrix
C = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        # also, C[i,j] = (1/m) * np.dot(M[:,i], M[:,j])
        C[i,j] = (1/m) * np.sum(M[:,i] * M[:,j])

# computing the eigenvalues and eigenvectors of this matrix
# since the covariance matrix is symmetric, the algorithm for Hermitian / symmetric
# matrices can be used
print C
values, vectors = np.linalg.eigh(C)
# TODO the vectors are the columns, right? keep that in mind
print values, vectors

print np.diag(C_ii)
# TODO what are the "true" eigenvalues / vectors? analytic formula? or just this?
true_values, true_vectors = np.linalg.eigh(np.diag(C_ii))

# using Euclidean distance to compare eigenvalues and eigenvectors
# TODO maybe angle would be more appropriate
# since scale of eigenvectors shouldn't matter
# TODO any scalar multiple of an eigenvector is also an eigenvector w/ same value, right?
# TODO sort first?
for v, v_hat in zip(true_values, values):
    print v, v_hat, (v - v_hat)**2

# TODO iterate over columns (assuming those are eigenvectors)
# TODO account for fact that vector could be pointing in either direction?
# sign of eigenvector doesnt matter, right?
# limit angle to be less than pi?
for v, v_hat in zip(true_vectors, vectors):
    radians = np.arccos(np.dot(v, v_hat) / \
        (np.linalg.norm(v) * np.linalg.norm(v_hat)))

    # TODO fix
    '''
    if radians > np.pi:
        radians = min(radians, 
    '''

    print v, v_hat, radians

plt.plot(M[:,0], M[:,1], '.')
# TODO equal axis limits to illustrate different ranges?

#plt.show()


