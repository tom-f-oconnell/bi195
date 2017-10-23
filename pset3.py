#!/usr/bin/env python

# TODO make a vim function to insert all this scientific python boilerplate
import numpy as np

from sympy import symbols
from sympy.plotting import plot
# the imaginary unit
from sympy import I
# seems to also work correctly on complex numbers (Euclidean norm of 2-vector)
from sympy.functions import Abs

import matplotlib.pyplot as plt

R, L, C, w = symbols('R L C w')

# TODO i might have thought this was going to be non-negative, but it seems
# that is not the case, unless I made a mistake?
hw_abs_squared = (w**2) / Abs(I*w + R/L - w**2*R*C)**2
hw_abs_squared_R = hw_abs_squared.subs(L, 1)
hw_abs_squared_R = hw_abs_squared_R.subs(C, 1)

# TODO log / log Bode plots?
plot(hw_abs_squared_R.subs(R, 1), (w, 0, 10), title='$R=1$')
plot(hw_abs_squared_R.subs(R, 5), (w, 0, 10), title='$R=5$')

plot(hw_abs_squared.subs(R, 1).subs(L, 5).subs(C, 1), (w, 0, 10), title='$L=5, C=1$')
plot(hw_abs_squared.subs(R, 1).subs(L, 1).subs(C, 5), (w, 0, 10), title='$L=1, C=5$')


