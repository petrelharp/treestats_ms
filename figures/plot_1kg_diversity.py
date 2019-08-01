import tskit, tszip, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

window_width = 5e5

popstyles = {"EAS" : (0, ()),
             "EUR" : (0, (1, 2)),
             "AFR" : (0, (3, 1)),
             "AMR" : (0, (3, 2, 1, 2)),
             "SAS" : (0, (3, 1, 3, 1, 1, 1))}
colors = plt.get_cmap("tab20").colors
pop_names = ['CHB', 'JPT', 'CHS', 'CDX', 'KHV', 'CEU', 'TSI', 'FIN', 'GBR', 'IBS', 'YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB', 'MXL', 'PUR', 'CLM', 'PEL', 'GIH', 'PJL', 'BEB', 'STU', 'ITU']
popcolors = { k : colors[j % len(colors)] for j, k in enumerate(pop_names) }

for chrom in range(1, 2):
    ts = tskit.load("1kg/1kg_chr{}.trees".format(chrom))

    pop_metadata = [json.loads(pop.metadata) for pop in ts.populations()]
    pn = [x['name'] for x in pop_metadata]
    for a,b in zip(pn, pop_names):
        assert(a == b)

    superpops = [x['super_population'] for x in pop_metadata]
    pop = [ts.node(ind.nodes[0]).population for ind in ts.individuals()]

    pop_nodes = [[] for _ in range(ts.num_populations)]
    for ind in ts.individuals():
        pop_nodes[pop[ind.id]].extend(ind.nodes)

    num_windows = int(ts.sequence_length / window_width)
    windows = np.linspace(0, ts.sequence_length, num_windows)
    x = windows[:-1] + np.diff(windows)/2

    bdiv = ts.diversity(pop_nodes, windows=windows, mode="branch")

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    for j in range(bdiv.shape[1]):
        dl = ax.plot(
                x, bdiv[:, j], label=pop_names[j],
                linestyle=popstyles[superpops[j]],
                color=popcolors[pop_names[j]])

    leg = ax.legend(
             fontsize = "small",
             loc = "upper left",
             handlelength=3,
             bbox_to_anchor = (1.01, 1.01),
             frameon = False,
             borderpad=0)
    fig.savefig("1kg/1kg_chr{}.diversity.pdf".format(chrom), bbox_inches = "tight")
    plt.close(fig)

