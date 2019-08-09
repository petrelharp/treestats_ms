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
statfile = f"{outbase}.f4.tsv"
plotfile = f"{outbase}.f4.pdf"

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


ts = tskit.load(treefile)

if ts.num_populations > 0:
    pop_metadata = [json.loads(pop.metadata) for pop in ts.populations()]
    for a, b, c in zip(pop_metadata, pop_names, superpops):
        assert(a['name'] == b)
        assert(a['super_population'] == c)

    pop = [ts.node(ind.nodes[0]).population for ind in ts.individuals()]

else:
    # this should be a Relate tree sequence
    assert(treefile.find("relate") >= 0)
    # in which case the metadata are in this external file, in order,
    # with one line per diploid
    pop = []
    with open("1kg/1000GP_Phase3_sub.sample", "r") as metafile:
        header = metafile.readline()
        assert(header == "ID POP GROUP SEX\n")
        for line in metafile:
            md = line.strip().split()
            for _ in range(2): # diploids!
                pop.append(pop_names.index(md[1]))

    
pop_nodes = [[] for _ in range(num_pops)]
for ind in ts.individuals():
    pop_nodes[pop[ind.id]].extend(ind.nodes)

# define the big windows
num_windows = int(ts.sequence_length / window_width)
windows = np.linspace(0, ts.sequence_length, num_windows + 1).astype('int')

try:
    with open(statfile, "r") as f:
        header = f.readline().strip().split()[1:]
    assert(len(header) == len(f4_pops))
    for a, b in zip(header, f4_pops):
        assert(a == ".".join(b))
    stat = np.loadtxt(statfile, skiprows=1)

except:
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

    indexes = [tuple([pop_names.index(u) for u in a]) for a in f4_pops]

    refined_stat = ts.f4(pop_nodes, indexes=indexes,
                         windows=refined_windows, mode=mode, span_normalise=False)
    # combine the many little windows back to the big ones
    stat = np.zeros((len(windows) - 1, len(indexes)))
    for k in range(len(indexes)):
        stat[:, k] = np.bincount(
                        refined_window_num[good_window],
                        weights=refined_stat[good_window, k],
                        minlength=num_windows)
        stat[nonmissing < min_window_content, k] = np.nan
        stat[nonmissing >= min_window_content, k] /= nonmissing[nonmissing >= min_window_content]

    np.savetxt(statfile, stat, 
               header="\t".join(".".join(a) for a in f4_pops),
               delimiter="\t")


#############
# plot results

x = windows[:-1] + np.diff(windows)/2

fig = plt.figure(figsize=(12 * len(windows) / 64, 6))
ax = fig.add_subplot(111)
for j in range(stat.shape[1]):
    dl = ax.plot(
            x, stat[:, j], label=stat_names[j],
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
