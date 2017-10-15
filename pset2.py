#!/usr/bin/env python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('dark')

def a(arraylike):
    return np.array(arraylike)


###############################################################################
# Problem 2
###############################################################################

print 'problem 2 (just for checking):'
A2 = a([[1,0,0],[0,3,6],[0,-1,-2]])

# each time i run the calculation, i come out with lambda = 0
# having algebraic multiplicity 1, and lambda = 1 having algebraic
# multiplicity 2
# YET the dimension of the null space of A2 - (1)*I is 1
# and it is 2 for (A2 (- 0)). (unless miscalculated). WHY?

def is_colinear(v1, v2):
    if np.isclose(np.dot(v1,v1)*np.dot(v2,v2),
        np.dot(v1,v2)**2):

        print '{} and {} were colinear'.format(v1, v2)
        return True

    else:
        print '{} and {} were NOT colinear'.format(v1, v2)
        return False


def check_eigenvalue(M, vec, val=None):
    vec = a(vec)
    print vec
    # TODO ...and why does only this eigenvector act like it 
    t = np.dot(M, vec)
    #assert np.allclose(t1, v1 * e1_manual), '\n{}\n{}'.format(v1, t1)
    print np.linalg.norm(t) / np.linalg.norm(vec)
    print t

#print 'eigenvalue for v1 (came from lambda=0 calculation): {}'.format(\

# checking A2 is idempotent numerically
assert np.allclose(A2, np.dot(A2, A2))

eivals, eivecs = np.linalg.eig(A2)
print 'eigenvalues and vectors determined computationally with ' + \
    'numpy:\neigenvalues:\n{}\neigenvectors:\n{}'.format(eivals, eivecs)

print 'matrix rank (perhaps determined through a similar ' + \
    'computation to above): {}'.format(np.linalg.matrix_rank(A2))

# eigenvectors i have determined manually
vs0 = {(1, 0, 0), (0, 1, -2)}
val0 = 0
for v in vs0:
    check_eigenvalue(A2, v, val=val0)

vs1 = {(0, 1, -3)}
val1 = 1
for v in vs1:
    check_eigenvalue(A2, v, val=val1)

# TODO use these to check
# from wolfram
correct_eivecs_wolfram = {(0,-3,1), (1,0,0), (0,-2,1)}
correct_eivecs_numpy = {tuple(eivecs[:,i]) for i in range(eivecs.shape[1])}
my_eivecs = vs0 | vs1
print correct_eivecs_wolfram
print correct_eivecs_numpy

for v1 in correct_eivecs_wolfram:
    assert sum(map(lambda x: is_colinear(v1, x), correct_eivecs_numpy)) == 1

for v1 in my_eivecs:
    assert sum(map(lambda x: is_colinear(v1, x), correct_eivecs_numpy)) == 1

###############################################################################
# Problem 3
###############################################################################

print '\nProblem 3 (just for checking):'
# j suffix creates an imaginary number
A3 = a([[0, -1j], [1j, 0]])
# would need sympy or similar to define x
# x3 = a([1, 'a'])
#eivals, eivecs = np.linalg.eig(A3)
#print 'eigenvalues and vectors determined computationally with ' + \
#    'numpy:\neigenvalues: {}\neigenvectors: {}'.format(eivals, eivecs)

# determined by hand
# for lambda = 1
eivec1 = a([1, -1j])
# for lambda = -1
eivec2 = a([1, 1j])

assert np.allclose(eivec1, (-1) * np.dot(A3, eivec1))
assert np.allclose(eivec2, np.dot(A3, eivec2))

# A = PDP-1

# eigenbasis
P = a([eivec1, eivec2]).T
Pinv = np.linalg.inv(P)

# could also just create D from the eigenvalues i calculated
D = np.dot(Pinv, np.dot(A3, P))

# change of basis to diagonalize A (and x)
# TODO come back to this
B = np.dot(P, np.dot(D, Pinv))
# TODO check rank is preserved? / other criteria for cob met?

def check_diagonal(M):
    # subtract M's diagonal from M, and check all zero
    assert np.allclose(M - np.diag(np.diag(M)), 0.0)

A3_diag = np.dot(B, A3)
check_diagonal(A3_diag)

print 'P:\n{}'.format(P)
# ?
print 'D:\n{}'.format(D)
print 'Pinv:\n{}\n'.format(Pinv)
print 'change of basis to diagonalize A (PD(P^-1)):\n{}\n'.format(B)
print 'A diagonalized:\n{}\n'.format(A3_diag)

# ?
D2 = a([[1, 0],[0, -1]])
A3_diag2 = np.dot(np.dot(P, np.dot(D2, Pinv)), A3)
print 'A diagonalized w/ D2:\n{}\n'.format(A3_diag2)


###############################################################################
# Problem 5
###############################################################################

def get_A5(x):
    return a([[x, -2], [3, 4]])

start = -15
stop = 30
num = 1000
atol_for_zero = 0.05

eigenvalues = []
xs = np.linspace(start, stop, num=num)
was_real_last = None
to_print = ''
for x in xs:
    eivals, _ = np.linalg.eig(get_A5(x))
    eigenvalues.append(eivals)

    # TODO tune so only one value flags w/ atol
    if np.isclose(0.0, eivals[0], atol=atol_for_zero):
        to_print += 'Eigenvalue 1 becomes zero at x={:.1f}\n'.format(x)
    if np.isclose(0.0, eivals[1], atol=atol_for_zero):
        to_print += 'Eigenvalue 2 becomes zero at x={:.1f}\n'.format(x)

    if np.any(np.nonzero(np.imag(eivals))):
        complex_evs = True
    else:
        complex_evs = False

    if was_real_last is None:
        x_on_last_change = x

    elif was_real_last and complex_evs:
        print 'Eigenvalues real on [{:.1f}, {:.1f}]'.format(x_on_last_change, x)
        x_on_last_change = x

    elif not (was_real_last or complex_evs):
        print 'Eigenvalues complex on [{:.1f}, {:.1f}]'.format(x_on_last_change, x)
        x_on_last_change = x

    was_real_last = False if complex_evs else True

if was_real_last:
    print 'Eigenvalues real on [{:.1f}, {:.1f}]'.format(x_on_last_change, x)
else:
    print 'Eigenvalues complex on [{:.1f}, {:.1f}]'.format(x_on_last_change, x)

print '\n{}'.format(to_print)

eigenvalues = np.stack(eigenvalues, axis=0)

fig = plt.figure()

def plot_complex_series(series, c, text):
    plt.plot(xs, np.real(series), c=c, label=('$Re$(' + text + ')'))
    plt.plot(xs, np.imag(series), c=c, linestyle='--', label=('$Im$(' + text + ')'))

plot_complex_series(eigenvalues[:,0], 'b', '$\lambda_1$')
plot_complex_series(eigenvalues[:,1], 'r', '$\lambda_2$')
plt.legend()
plt.title('Eigenvalues over parameter range')
plt.xlabel('$a$')
plt.ylabel('Eigenvalue component magnitude')

#plt.show()
