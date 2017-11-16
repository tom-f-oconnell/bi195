#!/usr/bin/env python

# Bi195 - Mathematics in Biology
# Problem set 6
# Thomas O'Connell

from __future__ import division
import numpy as np
import random
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('dark')

# Question 1: Polya's Urn

def urn_experiment(iteration_fn, initial_white, initial_black, \
    draws=10000, replicates=9):
    """
    """
    # can change to draws~=100 to better see effects
    # at earlier times
    n = int(np.ceil(np.sqrt(replicates)))
    fig, ax = plt.subplots(nrows=n, ncols=n, sharex=True, sharey=True)
    ax = ax.flatten()

    for j in range(replicates):
        num_white = initial_white
        num_black = initial_black
        fraction_black = np.zeros(draws) * np.nan

        for i in range(draws):
            dw, db = iteration_fn(num_white, num_black)
            num_white += dw
            num_black += db
            fraction_black[i] = num_black / (num_black + num_white)

        # TODO log x-axis?
        ax[j].plot(fraction_black, c='black')

    plt.suptitle('Fraction of black balls over time, {:,} iterations'.format(draws))
    fig.text(0.5, 0.04, 'Iteration', ha='center')
    fig.text(0.04, 0.5, 'Fraction of black balls', va='center', rotation='vertical')


# A
# add two of the same for each drawn
# start with one white and one black

def iteration_a(white, black):
    total = white + black
    # 0=white, 1=black
    ball = np.random.choice(2, p=[white/total, black/total])

    # effectively adding one ball of the same color
    if ball == 0:
        dw = 1
        db = 0
    elif ball == 1:
        dw = 0
        db = 1

    return dw, db

urn_experiment(iteration_a, 1, 1)

'''
There is more noise early in the process, because there are less total balls, 
so it is easier to change the proportion.

At late times, the fraction seems to remain steady around a value set early
in the process.

Each simulation does not lead to the same result. it seems possible to 
reach some kind of plateau at just about any value on [0,1].
'''


# B
# a model for genetic drift in populations of fixed size

# TODO try some w/ N = 20
N = 10

def iteration_b(white, black):
    total = white + black
    assert total == 2*N, 'total population not constant'

    # 0=white, 1=black
    ball1 = np.random.choice(2, p=[white/total, black/total])

    # to sample second ball without replacement
    if ball1 == 0:
        white -= 1
    else:
        black -= 1
    total -= 1
    ball2 = np.random.choice(2, p=[white/total, black/total])

    # if they are the same color, they are put back in
    if ball1 == ball2:
        dw = 0
        db = 0

    # otherwise, randomly pick one of the balls,
    # and put two of that color back in
    else:
        if np.random.choice([True, False]):
            dw = 1
            db = -1
        else:
            dw = -1
            db = 1

    return dw, db

urn_experiment(iteration_b, N, N, draws=1000)

'''
Unlike the previous experiment, here the total number of balls stays constant.
This means that the effect of any single draw does not diminish as the process
is repeated, so the fraction of any one color is similarly noisy throughout.

From replicates of this experiment, it seems inevitable that one color will
go extinct. In this model, there is no mechanism to favor the color with a 
lower fraction recovering from this state. It may eventually be unlikely to
draw one of the rare color, but when one is drawn, it is 50/50 as to whether
it will get subtracted from the population.
'''

plt.show()
