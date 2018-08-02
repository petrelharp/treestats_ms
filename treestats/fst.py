# Comparison of two methods for computing fst,
# one with genotypes, one with trees.
#
# These assume an infinite-sites model.

import msprime
import numpy as np


def variants_fst(tree_sequence, pop0, pop1, windows):
    '''
    Compute mean Fst across all variants in each window.
    '''
    n0 = len(pop0)
    n1 = len(pop1)
    num_windows = len(windows) - 1
    if windows[0] != 0.0:
        raise ValueError(
            "Windows must start at the start of the sequence (at 0.0).")
    if windows[-1] != tree_sequence.sequence_length:
        raise ValueError("Windows must extend to the end of the sequence.")
    for k in range(num_windows):
        if windows[k + 1] <= windows[k]:
            raise ValueError("Windows must be increasing.")
    S = np.zeros(num_windows)
    # index of *left-hand* end of the current window
    window_num = 0
    num_sites = 0
    for v in tree_sequence.variants():
        while v.position >= windows[window_num + 1]:
            if num_sites > 0:
                S[window_num] /= num_sites
            window_num += 1
            num_sites = 0
        num_sites += 1
        x0 = sum(v.genotypes[pop0])
        x1 = sum(v.genotypes[pop1])
        if (x0 + x1 > 0.0) and (x0 + x1 < n0 + n1):
            p0 = x0 / n0
            p1 = x1 / n1
            fst = 1 - (p0 * (1 - p0) + p1 * (1 - p1)) / ((p0 + p1) * (1 - (p0 + p1) / 2))
            S[window_num] += fst
    # do the last window
    if num_sites > 0:
        S[window_num] /= num_sites
    return S

def mutations_fst(tree_sequence, pop0, pop1, windows):
    '''
    Compute mean Fst across all variants in each window.
    Assumes an infinite-sites model.
    '''
    n0 = len(pop0)
    n1 = len(pop1)
    num_windows = len(windows) - 1
    if windows[0] != 0.0:
        raise ValueError(
            "Windows must start at the start of the sequence (at 0.0).")
    if windows[-1] != tree_sequence.sequence_length:
        raise ValueError("Windows must extend to the end of the sequence.")
    for k in range(num_windows):
        if windows[k + 1] <= windows[k]:
            raise ValueError("Windows must be increasing.")
    S = np.zeros(num_windows)
    # index of *left-hand* end of the current window
    window_num = 0
    num_sites = 0
    tt0 = tree_sequence.trees(tracked_samples=pop0,
                                   sample_counts=True)
    tt1 = tree_sequence.trees(tracked_samples=pop1,
                                   sample_counts=True)
    for t0, t1 in zip(tt0, tt1):
        for s in t0.sites():
            while s.position >= windows[window_num + 1]:
                if num_sites > 0:
                    S[window_num] /= num_sites
                window_num += 1
                num_sites = 0
            num_sites += 1
            for m in s.mutations:
                x0 = t0.num_tracked_samples(m.node)
                x1 = t1.num_tracked_samples(m.node)
                p0 = x0 / n0
                p1 = x1 / n1
                if (x0 + x1 > 0) and (x0 + x1 < n0 + n1):
                    fst = 1 - (p0 * (1 - p0) + p1 * (1 - p1)) / ((p0 + p1) * (1 - (p0 + p1) / 2))
                    S[window_num] += fst
    # do the last window
    if num_sites > 0:
        S[window_num] /= num_sites
    return S

