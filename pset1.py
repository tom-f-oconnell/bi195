#!/usr/bin/env python

from __future__ import division
import numpy as np

def differentiation_operator(n, dx, boundary_handling='none'):
    # TODO make this only calculate terms that would be fully
    # defined?
    A = np.zeros((n, n))
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            # assumes function starts at 0
            if i == (j - 1):
                A[i,j] = dx
            elif i == (j + 1):
                A[i,j] = -dx

    if boundary_handling == 'copy':
        if A.shape[1] >= 3:
            A[0,:] = A[1,:]
            A[-1,:] = A[-2,:]

    return A


# 5
def differentiate(f, dx=1, boundary_handling='copy', verbose=False):
    f = np.squeeze(f)
    # end of range (which starts at and includes 0. includes a.)
    a = f.shape[0] * dx
    # number of steps (# points = # steps + 1)
    n = f.size - 1

    A = differentiation_operator(f.size, n / (2 * a), boundary_handling=boundary_handling)
    g = np.dot(A, f)

    '''
    if boundary_handling == 'copy':
        # handle boundaries
        g[0] = g[1]
        g[-1] = g[-2]
    '''

    # 'correct' is what numpy calls this method of boundary handling
    if boundary_handling == 'correct':
        g = g[1:n]

    elif boundary_handling == 'none' or boundary_handling == 'copy':
        pass

    if verbose:
        print 'f', f
        print 'A:\n', A
        print 'g', g
        print ''

    return g


a = 4
n = 4
# testing with f(x)=x**2
f = np.linspace(0, a, num=(n+1)) ** 2
g = differentiate(f, verbose=True)

# is A invertible (5d)?
f1 = np.ones(5)
#f1[0] = 0
#f1[-1] = 0
#f1[3] = 1
c = 5
f2 = np.array(f1) * c
# TODO i thought this would work?
#f2[1:-1] = f2[1:-1] + c
#f2[3] = f2[3] + c
'''
g1 = differentiate(f1, boundary_handling='none', verbose=True)
g2 = differentiate(f2, boundary_handling='none', verbose=True)
'''
g1 = differentiate(f1, boundary_handling='copy', verbose=True)
g2 = differentiate(f2, boundary_handling='copy', verbose=True)

print 'f1', f1
print 'f2', f2
print 'g1', g1
print 'g2', g2
print ''

if np.allclose(g1, g2):
    print 'found a counterexample. ' + \
        'differentiation is not invertible!'
else:
    print 'differentiation may be invertible.'

A = differentiation_operator(2, 1, boundary_handling='copy')
print A

def det2by2(M):
    # should evaluate to 1 on I
    return M[0,0] * M[1,1] - M[0,1] * M[1,0]

print det2by2(np.eye(A.shape[0]))
print det2by2(A)
