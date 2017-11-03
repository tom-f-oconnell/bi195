#!/usr/bin/env python

from __future__ import division
import numpy as np
import re

def coinflip_sequences(length):
    if length == 0:
        return ['']
    else:
        sequences_one_shorter = coinflip_sequences(length - 1)
        return list(map(lambda s: s + 'H', sequences_one_shorter)) + \
            list(map(lambda s: s + 'T', sequences_one_shorter))

def analyze_sequences(seqs):
    hh_counts = np.array(map(lambda s: len(re.findall('(?=HH)', s)), seqs))
    #ht_counts = map(lambda s: len(re.findall('(?=HT)', s)), seqs)
    h_first_three = np.array(map(lambda s: s[:-1].count('H'), seqs))

    # using terminology from the (flawed) NYT article
    hh_percentage = hh_counts / h_first_three
    # gets the result in the NYT article
    print 'Mean HH-percentage, NYT method:', np.nanmean(hh_percentage)

    print 'Mean HH-percentage, considering # of H occurences:', \
        np.nansum(h_first_three * hh_percentage) / np.nansum(h_first_three)


for n in [4, 8, 12, 16]:
    analyze_sequences(coinflip_sequences(n))
    print ''
