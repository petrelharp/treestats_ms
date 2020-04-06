#!/usr/bin/env python3
import tskit, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys

usage = """
Usage:
    {} treefile maskfile mode window_width
Here `maskfile` should be a fasta of 0s and 1s.
""".format(sys.argv[0])

if len(sys.argv) != 5:
    raise ValueError(usage)

treefile = sys.argv[1]
maskfile = sys.argv[2]
mode = sys.argv[3]
window_width = float(sys.argv[4])

# min amount of data allowed to still deal with the window
min_window_content = window_width / 3

outbase = ".".join(treefile.split(".")[:-1]) + f".{mode}.{int(window_width)}"
divfile = f"{outbase}.diversity.tsv"
plotfile = f"{outbase}.diversity.pdf"
regionplotfile = f"{outbase}.region.diversity.pdf"


try:
    bdiv = np.loadtxt(divfile, skiprows=1)
    num_windows = bdiv.shape[0]
    windows = np.linspace(0, 63025522, num_windows + 1).astype('int')

except:

    ts = tskit.load(treefile)

    # define the big windows
    num_windows = int(ts.sequence_length / window_width)
    windows = np.linspace(0, ts.sequence_length, num_windows + 1).astype('int')

    pop_nodes = [[]]
    for ind in ts.individuals():
        pop_nodes[0].extend(ind.nodes)

    #############
    # get the mask

    mask = np.genfromtxt(maskfile, delimiter=1, skip_footer=1)
    masktail = np.genfromtxt(maskfile, delimiter=1, skip_header=len(mask))
    num_added_zeros = max(0, int(ts.sequence_length) - np.prod(mask.shape) - len(masktail))
    if num_added_zeros > 0:
        print(f"Adding {num_added_zeros} masked bases to the end of the mask.")

    mask = np.concatenate(
            [mask.reshape((mask.shape[0] * mask.shape[1],)),
             masktail,
             np.repeat(0, num_added_zeros)])

    mask = mask[:int(ts.sequence_length)]

    starts = np.where(np.diff(mask) > 0)[0]
    ends = np.where(np.diff(mask) < 0)[0]
    if len(ends) == len(starts) - 1:
        ends = np.concatenate([ends, np.array([len(mask) - 1])])
    lengths = ends - starts
    assert(np.all(lengths > 0))

    #############
    # compute the statistics in lots of little windows

    # which bases are in which window
    window_num = np.repeat(np.arange(num_windows), np.diff(windows))
    # find amount of each window that is not missing
    nonmissing = np.bincount(window_num, weights=mask)

    # find statistics in windows refined with segment breaks
    refined_windows = np.unique(np.concatenate([windows, starts, ends]))
    good_window = (mask[refined_windows[:-1]] == 1)
    refined_window_num = window_num[refined_windows[:-1]]

    refined_bdiv = ts.diversity(pop_nodes, windows=refined_windows,
                                mode=mode, span_normalise=False)

    # combine the many little windows back to the big ones
    bdiv = np.zeros((len(windows) - 1, len(pop_nodes)))
    for k in range(len(pop_nodes)):
        bdiv[:, k] = np.bincount(
                        refined_window_num[good_window],
                        weights=refined_bdiv[good_window, k],
                        minlength=num_windows)
        bdiv[nonmissing < min_window_content, k] = np.nan
        bdiv[nonmissing >= min_window_content, k] /= nonmissing[nonmissing >= min_window_content]


    # unmasked_bdiv = ts.diversity(pop_nodes, windows=windows, mode=mode)

    np.savetxt(divfile, bdiv, 
               header="diversity", delimiter="\t")

#############
# plot results

x = windows[:-1] + np.diff(windows)/2

fig = plt.figure(figsize=(12 * len(windows) / 64, 6))
ax = fig.add_subplot(111)
dl = ax.plot(x, bdiv)

fig.savefig(plotfile, bbox_inches = "tight")
plt.close(fig)
