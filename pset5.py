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

def evaluate_pc_estimate(C_ii, m, plot=False, verbose=False):
    """
    Args:
        m (int): number of samples / measurements
        C_ii (array-like): diagonal of covariance matrix.

        since covariance is diagonal, the data is restricted to have 
        independent components.

        plot (bool): if true, plots w/ sample data and estimates will be displayed.
        verbose (bool): if true, print extra stuff for debugging.

    Returns:
        value_error (float): sum of Euclidean distances between closest eigenvalues
        TODO check actually using closest eigenvalues

        vector_error (float): sum of angles between eigenvectors with closest
        eigenvalues. error angles restricted to [0, pi/2], since sign of
        eigenvectors should not matter. chose angle since normalization should
        also not matter.

        TODO should it be closest eigenvector or eigenvectors corresponding to closest
        eigenvalues? probably the latter? (what it actually is now, i suppose?)
    """
    # covariance diagonal
    # since these covariance matrices are zero off the diagonal
    # the r.v.s of each component are independent
    # TODO is that true? we could still have marginal normals and zero covariances,
    # but a non-normal joint distribution where the r.v.s are dependent, right?
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
    values, vectors = np.linalg.eigh(C)
    # TODO the vectors are the columns, right? keep that in mind
    if verbose:
        print C
        print values, vectors
        print np.diag(C_ii)

    # TODO what are the "true" eigenvalues / vectors? analytic formula? or just this?
    true_values, true_vectors = np.linalg.eigh(np.diag(C_ii))

    # using Euclidean distance to compare eigenvalues and eigenvectors
    # TODO maybe angle would be more appropriate
    # since scale of eigenvectors shouldn't matter
    # TODO any scalar multiple of an eigenvector is also an eigenvector w/ same value, right?
    # TODO sort first?
    value_error = 0
    for v, v_hat in zip(true_values, values):
        error = (v - v_hat)**2
        if verbose:
            print v, v_hat, error
        value_error += error

    # TODO iterate over columns (assuming those are eigenvectors)
    # TODO account for fact that vector could be pointing in either direction?
    # sign of eigenvector doesnt matter, right?
    # limit angle to be less than pi?
    vector_error = 0
    for v, v_hat in zip(true_vectors, vectors):
        # TODO see notes. will this handle pi + epsilon correctly? (ever return that?)
        radians = np.arccos(np.dot(v, v_hat) / \
            (np.linalg.norm(v) * np.linalg.norm(v_hat)))

        # TODO fix / check
        if radians > (np.pi / 2.0):
            radians = np.pi - radians

        vector_error += radians
        if verbose:
            print v, v_hat, radians

    if plot:
        plt.plot(M[:,0], M[:,1], '.')

    # TODO overlay actual and computed eigenvectors / values?
    # + illustration of similarity score?

    # TODO equal axis limits to illustrate different ranges?
    return value_error, vector_error


def evaluate_pc_estimates(C_ii, ms, nested_plot=False):
    """
    Args:
        C_ii (array-like): diagonal of covariance matrix. see above.
        ms (iterable): each element contains a number of samples. see above.

    Returns:
        value_errors (list): each element is a total eigenvalue error returned 
            for a given element of ms

        vector_errors (list): same, but for vector error as defined above.
    """
    value_errors = []
    vector_errors = []
    print 'C:\n', np.diag(C_ii)
    print 'm, value error, vector error'
    for m in ms:
        valerr, vecerr = evaluate_pc_estimate(C_ii, m, plot=nested_plot)
        value_errors.append(valerr)
        vector_errors.append(vecerr)

        # TODO truncate float printing
        print m, valerr, vecerr
    print ''

    # TODO do plotting / printing here?

    return value_errors, vector_errors

# 4.a
C_ii = [1, 0.25]
m = 1000
print '\n4a)'
evaluate_pc_estimate(C_ii, m, plot=True, verbose=True)

# 4.b
ms = [30, 100, 300, 1000]
print '\n4b)'
evaluate_pc_estimates(C_ii, ms)
# TODO plot?

# 4.c
C_ii = [1, 0.9]
print '4c)'
evaluate_pc_estimates(C_ii, ms)

# TODO same progression of m?
# 4.d
# TODO what is this supposed to expose? is this the "curse of dimensionality"?
C_ii = [1] + ([0.25] * 4)
print '4d)'
evaluate_pc_estimates(C_ii, ms)

# TODO statistically test whether error is significantly higher (for given m)
# or decreasing slower, at higher dimensions?

#plt.show()

