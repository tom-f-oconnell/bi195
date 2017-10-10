#!/usr/bin/env python

from __future__ import division
import numpy as np

'''
Problem 5
'''

def differentiation_operator(n, dx, boundary_handling='none'):
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


def differentiate(f, dx=1, boundary_handling='copy', verbose=False):
    f = np.squeeze(f)
    # end of range (which starts at and includes 0. includes a.)
    a = f.shape[0] * dx
    # number of steps (# points = # steps + 1)
    n = f.size - 1

    A = differentiation_operator(f.size, n / (2 * a), boundary_handling=boundary_handling)
    g = np.dot(A, f)

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

print 'Problem 5:'

a = 4
n = 4
# testing with f(x)=x**2
f = np.linspace(0, a, num=(n+1)) ** 2
g = differentiate(f, verbose=True)

# is A invertible (5d)?
f1 = np.ones(5)
c = 5
f2 = np.array(f1) * c
g1 = differentiate(f1, boundary_handling='copy', verbose=True)
g2 = differentiate(f2, boundary_handling='copy', verbose=True)

print 'f1', f1
print 'f2', f2
print 'g1', g1
print 'g2', g2

if np.allclose(g1, g2):
    print 'found a counterexample. ' + \
        'differentiation is not invertible!'
else:
    print 'differentiation may be invertible.'

A = differentiation_operator(2, 1, boundary_handling='copy')
print '\nthe differentiation operator under test:'
print A

def det2by2(M):
    return M[0,0] * M[1,1] - M[0,1] * M[1,0]

# should evaluate to 1 on I
assert np.isclose(det2by2(np.eye(A.shape[0])), 1.0)
print 'its determinant:'
print det2by2(A)


'''
Problem 6
'''

def leading_coefficient_idx(v):
    """
    Doesn't check whether v is a zero vector.
    """
    return np.argmax(np.abs(v) > 1e-6)


def gaussian_elimination(A, stop='gauss', det_change_factor=False):
    """
    Returns the row echelon form of A.

    Goal:
    -transform matrix through [row XOR col] operations to one that is in [row XOR col]
     echelon form

    Uses:
    -calculating determinants
    -"solving" systems of equations

    Operations available ("elementary row operations"):
    1: swap two rows
    2: multiply a row by a nonzero scalar
    3: add to one row a scalar multiple of another (seems to unecessarily include rule 2, 
       doesn't specifcy nonzero here either, on Wiki, but I think they mean to)

    -should run in O(n^3)
    TODO where n is number of elements or one dimension of (square?) matrix?
    -sometimes unstable, but generally not (for some classes of matrices at least)
    """
    if stop not in {'gauss', 'gauss-jordan'}:
        raise ValueError('invalid stopping condition')
    
    def enumerate_leads(M):
        return [(row_idx, leading_coefficient_idx(row)) for row_idx, row in enumerate(M)]

    def sort_rows_by_leading_index(M):
        pairs = enumerate_leads(M)
        # sort output places keys from small to big
        sorted_row_indices = [x for x, y in sorted(pairs, key=lambda x: x[1])]
        return M[sorted_row_indices, :]
    
    zero_row = np.zeros(A.shape[1])
    def is_zeros(row):
        return np.allclose(row, zero_row)

    if det_change_factor:
        d = 1

    A = sort_rows_by_leading_index(A)
    curr_row_idx = 0

    while curr_row_idx < A.shape[0]:
        curr_lead = leading_coefficient_idx(A[curr_row_idx])

        # find one row with same index for its leading coefficient
        # if none exist, continue to next row
        other_leads = enumerate_leads(A[curr_row_idx + 1:, :])
        ties = [i + 1 + curr_row_idx for i, l in other_leads if l == curr_lead]

        # make it so no leading coefficients have the same index
        while len(ties) != 0:
            i = ties.pop()
            scale = A[i, curr_lead] / A[curr_row_idx, curr_lead]
            A[i,:] = A[i,:] - scale * A[curr_row_idx, :]

        A[curr_row_idx + 1:, :] = sort_rows_by_leading_index(A[curr_row_idx + 1:, :])
        curr_row_idx += 1

    if stop == 'gauss-jordan':
        lead_indices = []
        # make all leading coefficients 1
        # by dividing the rows by themselves
        for i in range(A.shape[0]):
            row = A[i,:]
            # assumes sorted s.t. all zero rows beneath all nonzero rows
            if is_zeros(row):
                break
            lead_idx = leading_coefficient_idx(row)
            lead_indices.append((i,lead_idx))
            leading_coeff = A[i, lead_idx]
            A[i,:] = row / leading_coeff

        # make all other entries in columns with leading indices zero
        # by subtracting (multiples of) the row with the lead from the others
        for i, j in lead_indices:
            lead = A[i,j]
            for i_other in range(A.shape[0]):
                if i != i_other and not np.isclose(A[i_other,j], 0):
                    scale = A[i_other,j] / A[i,j]
                    A[i_other,:] = A[i_other,:] - scale * A[i,:]
    return A


def solve(A, b):
    """
    Returns the x from Ax=b
    """
    # make augmented
    Ab = np.concatenate((A, b), axis=1)
    S = gaussian_elimination(Ab, stop='gauss-jordan')
    x = S[:,-1]
    return x


print '\nProblem 6:'
A = np.array( \
    [[2, -1, -3], \
     [-2, 5, 1], \
     [1, 4, 6]])
b = np.array([-21.25, 41.5833, 35]).reshape(3, 1)

print 'A:\n', A
print 'b:\n', b
print 'x:\n', solve(A, b).reshape(3, 1)
print 'numpy.linalg.solve(A,b), for comparison:\n', \
    np.linalg.solve(A, b)
