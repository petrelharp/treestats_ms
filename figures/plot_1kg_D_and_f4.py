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

outbase = ".".join(treefile.split(".")[:-1]) + f".{mode}.{int(window_width)}"
statfile = f"{outbase}.divergences.tsv" # created by calc_1kg_divergences.py
D_plotfile = f"{outbase}.D.pdf"
f4_plotfile = f"{outbase}.f4.pdf"

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

# compute D, f4:
dstat = np.zeros((len(windows) - 1, len(f4_pops)))
f4stat = np.zeros((len(windows) - 1, len(f4_pops)))
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
        f4stat[:, m] = numer / nonmissing


#############
# plot results

x = windows[:-1] + np.diff(windows)/2

# plot D
fig = plt.figure(figsize=(12 * np.sqrt(len(windows) / 64), 6))
ax = fig.add_subplot(111)
for j in range(dstat.shape[1]):
    dl = ax.plot(
            x, dstat[:, j], label=stat_names[j],
            color=stat_colors[stat_names[j]])

leg = ax.legend(
         fontsize = "small",
         loc = "upper left",
         handlelength=3,
         bbox_to_anchor = (1.01, 1.01),
         frameon = False,
         borderpad=0)
fig.savefig(D_plotfile, bbox_inches = "tight")
plt.close(fig)


# and f4
fig = plt.figure(figsize=(12 * np.sqrt(len(windows) / 64), 6))
ax = fig.add_subplot(111)
for j in range(f4stat.shape[1]):
    dl = ax.plot(
            x, f4stat[:, j], label=stat_names[j],
            color=stat_colors[stat_names[j]])

leg = ax.legend(
         fontsize = "small",
         loc = "upper left",
         handlelength=3,
         bbox_to_anchor = (1.01, 1.01),
         frameon = False,
         borderpad=0)
fig.savefig(f4_plotfile, bbox_inches = "tight")
plt.close(fig)



