import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize


def tskit_model(n, a, b, c):
    # The model is based on the proof that we have O(n + rho log n) edges
    # and we look at each of these twice as we go through the trees.

    # TODO not quite sure about this, need to think about it again.
    # Excellent fit though.
    rho = 0.025  # Value from the file
    z = rho * np.log(n)
    # TODO: did KT mess this up?
    return a*n + b*(z) * np.log(n) + c
    # return a * (n+z) + b * np.log(n) + c


symbols = {"tskit": "o", "allel": "^", "libseq": "s"}

for infile in ["benchmarks_without_copy_longer_genome.txt"]:
    b = pd.read_csv(infile, sep=' ')
    b = b.round({'Nu': 4})
    b['seconds_snp'] = b.seconds/b.nmutations
    grp = b.groupby(['toolkit'])
    fig = plt.figure()
    for n, g in grp:
        if n == 'tskit':
            line, = plt.loglog(g.nsam, g.seconds_snp, symbols[n], label=n)
            fit_params, _ = optimize.curve_fit(
                tskit_model, g.nsam, g.seconds_snp)
            fit = tskit_model(g.nsam, *fit_params)
            plt.loglog(g.nsam, fit, color=line.get_color(),
                       label='__nolegend__')
        else:
            line, = plt.loglog(g.nsam, g.seconds_snp, symbols[n], label=n)
            plt.loglog(g.nsam, g.seconds_snp, linestyle=':', color=line.get_color(),
                       label='__nolegend__')
    plt.xlabel("Sample size (haploid)")
    plt.ylabel("Run time per variant (seconds)")
    of = infile.replace("txt", "pdf")
    plt.legend(loc='upper left', fontsize='x-small')
    plt.tight_layout()
    plt.savefig(of)
    plt.clf()
