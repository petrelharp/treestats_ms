#!/usr/bin/env python3
import tskit, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys

usage = """
Usage:
    {} treefile maskfile
Here `maskfile` should be a fasta of 0s and 1s.
""".format(sys.argv[0])


if len(sys.argv) != 3:
    raise ValueError(usage)

treefile = sys.argv[1]
maskfile = sys.argv[2]

outbase = ".".join(treefile.split(".")[:-1]) + f".node"
plotfile = f"{outbase}.f4.pdf"
timeplotfile = f"{outbase}.time.f4.pdf"

method_label = {
        "1kg/relate_chr20.trees" : "Relate",
        "1kg/tgp_geva_chr20.trees" : "GEVA",
        "1kg/1kg_chr20.trees" : "tsinfer"
        }[treefile]

pop_names = ['CHB', 'JPT', 'CHS', 'CDX', 'KHV', 'CEU', 'TSI', 'FIN', 'GBR', 'IBS', 'YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB', 'MXL', 'PUR', 'CLM', 'PEL', 'GIH', 'PJL', 'BEB', 'STU', 'ITU']
superpops = ['EAS', 'EAS', 'EAS', 'EAS', 'EAS', 'EUR', 'EUR', 'EUR', 'EUR', 'EUR', 'AFR', 'AFR', 'AFR', 'AFR', 'AFR', 'AFR', 'AFR', 'AMR', 'AMR', 'AMR', 'AMR', 'SAS', 'SAS', 'SAS', 'SAS', 'SAS']
num_pops = len(pop_names)

# which populations?
f4_pops = [('PUR', 'TSI', 'GWD', 'JPT'),
           ('ASW', 'CEU', 'MSL', 'CHB'),
           ('IBS', 'GBR', 'FIN', 'JPT'),
           ('TSI', 'CEU', 'FIN', 'CHB')]
f4_indexes = [tuple(pop_names.index(i) for i in a) for a in f4_pops]

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


# bf4 = ts.f4(pop_nodes, indexes=f4_indexes,
#             mode="branch")[0, :]

nf4 = ts.f4(pop_nodes, indexes=f4_indexes,
           mode="node", span_normalise=False)[0, :, :]

if method_label == "Relate":
    num_bins = 500
else:
    num_bins = 50

node_times = ts.tables.nodes.time
time_bins = np.concatenate([[0], np.exp(np.linspace(np.log(100), np.log(max(node_times) + 1), num_bins - 1))])
node_bins = np.digitize(node_times, time_bins)

f4_binned = np.zeros((num_bins, nf4.shape[1]))
for j in range(nf4.shape[1]):
    f4_binned[:, j] = np.bincount(node_bins, weights=nf4[:, j], minlength=num_bins)
    # f4_binned[:, j] = (np.bincount(node_bins, weights=nf4[:, j], minlength=num_bins)
    #                    / np.bincount(node_bins, minlength=num_bins))

f4_binned /= np.nansum(f4_binned, axis=0)


#### plot it

fig = plt.figure(figsize=(6, 3))
ax = fig.add_subplot(111)
ax.set_xlabel("time ago")
ax.set_ylabel("proportion of node f4")
# if method_label == "Relate":
#     ax.set_xlim(0, 20000)
# elif method_label == "GEVA":
#     ax.set_xlim(0, 10000)

for j, sn in enumerate(stat_names):
    fl = ax.semilogx(time_bins, f4_binned[:, j], label=sn)

ax.legend(
        fontsize = "small",
        frameon = False,
        borderpad = 0)

fig.savefig(plotfile, bbox_inches = "tight")
plt.close(fig)

###### load and plot branch-length-weighted f4 values

time_f4 = np.zeros((num_bins, len(f4_indexes)))
for j, f4p in enumerate(f4_pops):
    sn = ".".join(f4p)
    datafile = f"{outbase}.{sn}.realf4.binned.txt"
    rf4 = np.loadtxt(datafile)
    time_f4[:, j] = rf4 / np.nansum(rf4)


fig = plt.figure(figsize=(6, 3))
ax = fig.add_subplot(111)
ax.set_xlabel("time ago (generations)")
ax.set_ylabel("normalised branch f4")

for j, sn in enumerate(stat_names):
    fl = ax.semilogx(time_bins, time_f4[:, j], label=sn)

ax.legend(
        fontsize = "small",
        frameon = False,
        borderpad = 0)

fig.savefig(timeplotfile, bbox_inches = "tight")
plt.close(fig)

