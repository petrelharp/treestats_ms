import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize


def tskit_model(logn, a, b, c):
    n = np.exp(logn)
    rho = 0.025  # Value from the file
    return np.log(a * n + b * rho * logn * logn + c)


symbols = {"tskit": "o", "allel": "^", "libseq": "s"}
infile = "benchmarks_without_copy_longer_genome.txt"

b = pd.read_csv(infile, sep=' ')
b = b.round({'Nu': 4})
b['seconds_snp'] = b.seconds/b.nmutations
b['snp_rate'] = 1 / b.seconds_snp
grp = b.groupby(['toolkit'])
fig = plt.figure()
for n, g in grp:
    if n == 'tskit':
        line, = plt.loglog(g.nsam, g.seconds_snp, symbols[n], label=n)
        # We fit the model in log-space.
        fit_params, _ = optimize.curve_fit(
            tskit_model, np.log(g.nsam), np.log(g.seconds))
        fit = np.exp(tskit_model(np.log(g.nsam), *fit_params)) / g.nmutations
        plt.loglog(g.nsam, fit, color=line.get_color(),
                   label='__nolegend__')

        k = 12
        X = g.nsam.values[:]
        x = X[k] + (X[k + 1] - X[k]) / 2
        Y = fit.values[:]
        y = Y[k] + (Y[k + 1] - Y[k]) / 2
        plt.annotate(
            "Fitted model: $an + b\\rho\\log(n)^2 + c$",
            xy=(x, y),
            textcoords="offset pixels", xytext=(-85, -30),
            arrowprops={"arrowstyle": "->"})
    else:
        line, = plt.loglog(g.nsam, g.seconds_snp, symbols[n], label=n)
        plt.loglog(g.nsam, g.seconds_snp, linestyle=':', color=line.get_color(),
                   label='__nolegend__')
    yoffset = 0
    if n == "allel":
        yoffset = 5
    row = g.iloc[-1]
    rate = "{:,} var/s".format(int(row.snp_rate))
    plt.annotate(
        rate, xy=(row.nsam, row.seconds_snp), textcoords="offset pixels",
        xytext=(20, yoffset))

plt.xlabel("Sample size (haploid)")
plt.ylabel("Run time per variant (seconds)")
of = infile.replace("txt", "pdf")
plt.legend(loc='upper left')
plt.tight_layout(rect=(0, 0, 0.85, 1))
plt.savefig(of)
plt.clf()
