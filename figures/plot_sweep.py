import pyslim, tskit, msprime
import numpy as np
import os, sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

usage = """
Usage:
    {} treefile Ne
""".format(sys.argv[0])

if len(sys.argv) != 3:
    raise ValueError(usage)

treefile = sys.argv[1]
Ne = int(sys.argv[2])

seed = 23
mut_rate = 1e-9

outbase = ".".join(treefile.split(".")[:-1]) + ".{}".format(mut_rate)
newtrees = outbase + ".recap.trees"

try:
    ts = pyslim.load(newtrees)
except:
    ots = pyslim.load(treefile)
    rts = ots.recapitate(1e-8, Ne=Ne)
    mts = pyslim.SlimTreeSequence(msprime.mutate(rts, rate=mut_rate, random_seed=seed, keep=True))
    ts = mts.simplify()
    ts.dump(newtrees)

# get swept mutations and selection coefficients
swept_id = np.where(np.diff(ts.tables.mutations.metadata_offset) > 0)[0]
swept = [ts.mutation(m) for m in swept_id]
s = np.array([mut.metadata[0].selection_coeff for mut in swept])
spos = np.array([ts.site(mut.site).position for mut in swept])

windows = np.arange(0, ts.sequence_length, 5e5)
windows[-1] = ts.sequence_length
x = windows[:-1] + np.diff(windows)/2
site_div = ts.diversity([ts.samples()], windows=windows, mode='site')
branch_div = mut_rate * ts.diversity([ts.samples()], windows=windows, mode='branch')

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(111)
maxy = 1.1 * max(np.max(np.abs(site_div)), np.max(np.abs(branch_div)))
ax.set_ylim(0.0, maxy)
sl, = ax.plot(x, site_div,
             alpha=1.0,
             color='black',
             label="site")
bl, = ax.plot(x, branch_div,
             alpha=1.0,
             color='red',
             label="branch")
sl = ax.axvline(spos[0], linestyle=":", label="sweep")
for k, (pos, sval) in enumerate(zip(spos, s)):
    if k > 0:
        ax.axvline(pos, linestyle=":")
    # # text labels with selection coefficients?
    # at = ax.text(pos, 0.95*maxy, "s={:.2}".format(sval))

ax.set_xlabel("chromosome position (bp)")
ax.set_ylabel("mean diversity (pi)")
ax.legend()
fig.savefig(outbase + ".diversity.pdf", bbox_inches = "tight")
plt.close(fig)

