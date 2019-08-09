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
Here `maskfile` is unused.
""".format(sys.argv[0])


if len(sys.argv) != 5:
    raise ValueError(usage)

treefile = sys.argv[1]
maskfile = sys.argv[2]
mode = sys.argv[3]
window_width = float(sys.argv[4])

outbase = ".".join(treefile.split(".")[:-1]) + f".{mode}.{int(window_width)}"
statfile = f"{outbase}.divergences.tsv" # created by calc_1kg_divergences.py
D_plotfile = f"{outbase}.D.jackknife.pdf"
f4_plotfile = f"{outbase}.f4.jackknife.pdf"

pop_names = ['CHB', 'JPT', 'CHS', 'CDX', 'KHV', 'CEU', 'TSI', 'FIN', 'GBR', 'IBS', 'YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB', 'MXL', 'PUR', 'CLM', 'PEL', 'GIH', 'PJL', 'BEB', 'STU', 'ITU']
superpops = ['EAS', 'EAS', 'EAS', 'EAS', 'EAS', 'EUR', 'EUR', 'EUR', 'EUR', 'EUR', 'AFR', 'AFR', 'AFR', 'AFR', 'AFR', 'AFR', 'AFR', 'AMR', 'AMR', 'AMR', 'AMR', 'SAS', 'SAS', 'SAS', 'SAS', 'SAS']
num_pops = len(pop_names)

# which populations?
f4_pops = [('PUR', 'TSI', 'GWD', 'CHB'),
           ('ASW', 'CEU', 'MSL', 'CHB'),
           ('GBR', 'IBS', 'ITU', 'JPT')]

stat_names = [f"{a},{b};{c},{d}" for a, b, c, d in f4_pops]
colors = plt.get_cmap("tab20").colors
stat_colors = { k : colors[j] for j, k in enumerate(stat_names) }

# read in already-computed divergences
all_pairs = [(i, j) for i in range(num_pops) for j in range(i, num_pops)]
all_pairs_names = [tuple([pop_names[i] for i in a]) for a in all_pairs]

try:
    header = ["start", "end", "nonmissing"] + [".".join([pop_names[i] for i in a]) for a in all_pairs]
    with open(statfile, "r") as f:
        the_header = f.readline().strip().split()[1:]
    assert(len(header) == len(the_header))
    for a, b in zip(header, the_header):
        assert(a == b)
    divs = np.loadtxt(statfile, skiprows=1)
except:
    raise ValueError("Need to calculate divergences first, with calc_1kg_divergences.py.")


# get windows from start, end
num_windows = divs.shape[0]
assert(np.all(divs[1:,0] == divs[:-1,1]))
windows = np.concatenate([divs[:, 0], [divs[-1, 1]]])
nonmissing = divs[:, 2]
divs = divs[:, 3:]


def op(a,b):
    return (a,b) if pop_names.index(a) < pop_names.index(b) else (b,a)

def bincount_array(x, weights, minlength=0):
    if minlength == 0:
        minlength = max(x) + 1
    out = np.zeros((minlength, weights.shape[1]))
    for k in range(weights.shape[1]):
        out[:, k] = np.bincount(x, weights[:, k], minlength=minlength)
    return out

def calc_D(divs, nonmissing):
    dstat = np.zeros((divs.shape[0], len(f4_pops)))
    for m, (i, j, k, l) in enumerate(f4_pops):
        ab = all_pairs_names.index(op(i,j))
        ac = all_pairs_names.index(op(i,k))
        ad = all_pairs_names.index(op(i,l))
        bc = all_pairs_names.index(op(j,k))
        bd = all_pairs_names.index(op(j,l))
        cd = all_pairs_names.index(op(k,l))
        numer = (divs[:, ad] + divs[:, bc] - divs[:, ac] - divs[:, bd])
        denom = (2 * divs[:, ab] + 2 * divs[:, cd]
                 - divs[:, ad] - divs[:, bc] - divs[:, ac] - divs[:, bd])
        with np.errstate(divide='ignore', invalid='ignore'):
            dstat[:, m] = numer / denom
    return dstat

def calc_f4(divs, nonmissing):
    stat = np.zeros((divs.shape[0], len(f4_pops)))
    for m, (i, j, k, l) in enumerate(f4_pops):
        ab = all_pairs_names.index(op(i,j))
        ac = all_pairs_names.index(op(i,k))
        ad = all_pairs_names.index(op(i,l))
        bc = all_pairs_names.index(op(j,k))
        bd = all_pairs_names.index(op(j,l))
        cd = all_pairs_names.index(op(k,l))
        numer = (divs[:, ad] + divs[:, bc] - divs[:, ac] - divs[:, bd])
        with np.errstate(divide='ignore', invalid='ignore'):
            stat[:, m] = numer / nonmissing
    return stat

def jackknife(divs, fun, blocksize=10):
    """
    Calculate fun on divs after dropping each of the rows in each block.
    """
    num_blocks = int(np.ceil(divs.shape[0] / blocksize))
    bindex = np.repeat([np.arange(blocksize)], num_blocks, 0)
    bindex = bindex.reshape((np.prod(bindex.shape),))[:divs.shape[0]]
    block = np.repeat(np.arange(num_blocks), blocksize)[:divs.shape[0]]
    out = []
    for k in range(blocksize):
        subset = (bindex != k)
        sub_divs = bincount_array(block[subset], divs[subset, :], minlength=num_blocks)
        sub_nonmissing = np.bincount(block[subset], nonmissing[subset], minlength=num_blocks)
        out.append(fun(sub_divs, sub_nonmissing))
    return np.array(out / np.sqrt(blocksize))


num_stats = len(f4_pops)
blocksizes = np.arange(5, 60, 5)

##########
# D statistic

for statname, statfun, plotfile in [
     ["D", calc_D, D_plotfile],
     ["f4", calc_f4, f4_plotfile]]:

    means = np.zeros((len(blocksizes), num_stats))
    sds = np.zeros((len(blocksizes), num_stats))

    for j, blocksize in enumerate(blocksizes):
        jD = jackknife(divs, statfun, blocksize=blocksize)
        # average jackknife SD across windows
        sds[j, :] = np.nanmean(np.std(jD, axis=0) / blocksize, axis=0)
        # average jackknife estimate across windows
        means[j, :] = np.nanmean(np.mean(jD, axis=0) / blocksize, axis=0)

    big_windows = np.array([b * round(np.mean(np.diff(windows)), -2) for b in blocksizes])

    fig = plt.figure(figsize=(12, 6))

    ax = fig.add_subplot(121)
    ax.set_xlabel("window size")
    ax.set_ylabel(f"mean {statname} across windows")
    for j in range(num_stats):
        dl = ax.plot(
                big_windows, means[:, j], label=stat_names[j],
                color=stat_colors[stat_names[j]])

    ax = fig.add_subplot(122)
    ax.set_xlabel("window size")
    ax.set_ylabel(f"mean jackknife SD of {statname} across windows")
    for j in range(num_stats):
        dl = ax.plot(
                big_windows, sds[:, j], label=stat_names[j],
                color=stat_colors[stat_names[j]])

    leg = ax.legend(
             fontsize = "small",
             loc = "upper left",
             handlelength=3,
             bbox_to_anchor = (1.01, 1.01),
             frameon = False,
             borderpad=0)
    fig.savefig(plotfile, bbox_inches = "tight")
    plt.close(fig)
