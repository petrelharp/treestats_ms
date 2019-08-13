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



#########################

# real_nf4 = node_f4(ts, pop_nodes, indexes=f4_indexes,
#                    span_normalise=False)

def node_f4(ts, sample_sets, indexes, windows=None, span_normalise=True):
    windows = ts.parse_windows(windows)
    out = np.zeros((len(windows) - 1, ts.num_nodes, len(indexes)))
    node_times = ts.tables.nodes.time
    for i, (iA, iB, iC, iD) in enumerate(indexes):
        A = sample_sets[iA]
        B = sample_sets[iB]
        C = sample_sets[iC]
        D = sample_sets[iD]
        tA = len(A)
        tB = len(B)
        tC = len(C)
        tD = len(D)
        denom = np.float64(tA * tB * tC * tD)
        for j in range(len(windows) - 1):
            begin = windows[j]
            end = windows[j + 1]
            S = np.zeros(ts.num_nodes)
            for k, (t1, t2, t3, t4) in enumerate(zip(ts.trees(tracked_samples=A),
                                      ts.trees(tracked_samples=B),
                                      ts.trees(tracked_samples=C),
                                      ts.trees(tracked_samples=D))):
                print(k)
                if t1.interval[1] <= begin:
                    continue
                if t1.interval[0] >= end:
                    break
                SS = np.zeros(ts.num_nodes)
                for u in t1.nodes():
                    h = node_times[t1.parent(u)] - node_times[u]
                    # count number of pairwise paths going through u
                    nA = t1.num_tracked_samples(u)
                    nB = t2.num_tracked_samples(u)
                    nC = t3.num_tracked_samples(u)
                    nD = t4.num_tracked_samples(u)
                    # ac|bd - ad|bc
                    SS[u] += h * (nA * nC * (tB - nB) * (tD - nD)
                                  + (tA - nA) * (tC - nC) * nB * nD)
                    SS[u] -= h * (nA * nD * (tB - nB) * (tC - nC)
                                  + (tA - nA) * (tD - nD) * nB * nC)
                S += SS*(min(end, t1.interval[1]) - max(begin, t1.interval[0]))
            out[j, :, i] = S / denom
            if span_normalise:
                out[j, :, i] /= (end - begin)
    return out

for f4p, f4i in zip(stat_names, f4_indexes):
    s_all_nodes = []
    s_pop_nodes = []
    s_nm = {}
    s_pm = {}
    i = 0
    n = 0
    for j in f4i:
        s_pm[j] = i
        i += 1
        pn = pop_nodes[j]
        new_pn = []
        for u in pn:
            s_all_nodes.append(u)
            s_nm[u] = n
            new_pn.append(n)
            n += 1
        s_pop_nodes.append(new_pn)

    s_ts = ts.simplify(s_all_nodes)

    nf4 = node_f4(s_ts, s_pop_nodes, indexes=[(0, 1, 2, 3)],
                       span_normalise=False).reshape((s_ts.num_nodes,))

    if method_label == "Relate":
        num_bins = 500
    else:
        num_bins = 50

    node_times = s_ts.tables.nodes.time
    time_bins = np.concatenate([[0], np.exp(np.linspace(np.log(100), np.log(max(node_times) + 1), num_bins - 1))])
    node_bins = np.digitize(node_times, time_bins)
    f4_binned = np.bincount(node_bins, weights=nf4, minlength=num_bins)

    sn = ".".join(f4p)
    outfile = f"{outbase}.{sn}.realf4.txt"
    np.savetxt(outfile, np.column_stack([node_times, nf4]), header=f"time\t{sn}")
    sumfile = f"{outbase}.{sn}.realf4.binned.txt"
    np.savetxt(sumfile, f4_binned, header=sn)
