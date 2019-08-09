#!/usr/bin/env python3
import pyslim, tskit, msprime
import numpy as np
import os, sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

usage = """
Usage:
    {} treefile Ne mut_rate
""".format(sys.argv[0])

if len(sys.argv) != 4:
    raise ValueError(usage)

treefile = sys.argv[1]
Ne = int(sys.argv[2])
mut_rate = float(sys.argv[3])

nreps = 20

outbase = ".".join(treefile.split(".")[:-1])
newtrees = outbase + ".recap.trees"

try:
    ts = pyslim.load(newtrees)
except:
    ots = pyslim.load(treefile)
    rts = ots.recapitate(1e-8, Ne=Ne)
    ts = rts.simplify()
    ts.dump(newtrees)

# get swept mutations and selection coefficients
s = np.array([mut.metadata[0].selection_coeff for mut in ts.mutations()])
spos = np.array([ts.site(mut.site).position for mut in ts.mutations()])
freqs = np.array([ts.at(ts.site(mut.site).position).num_samples(mut.node) for mut in ts.mutations()])
spos = spos[freqs > 0.5 * ts.num_samples]
s = s[freqs > 0.5 * ts.num_samples]

windows = np.arange(0, ts.sequence_length, 5e5)
windows[-1] = ts.sequence_length
x = windows[:-1] + np.diff(windows)/2

########################
# Plot lines for 

all_branch_div = ts.diversity([ts.samples()], windows=windows, mode='branch')

# independent versions of the same rate, low

low_mut_rate = mut_rate
all_site_div = np.zeros((len(windows) - 1, nreps))
for k in range(nreps):
    mts = msprime.mutate(ts, rate=low_mut_rate, keep=True)
    all_site_div[:, k] = (1/low_mut_rate) * mts.diversity([ts.samples()], windows=windows, mode='site')[:, 0]

# and independent versions of the same rate, high

high_mut_rate = mut_rate * 10
all_site_div_high = np.zeros((len(windows) - 1, nreps))
for k in range(nreps):
    mts = msprime.mutate(ts, rate=high_mut_rate, keep=True)
    all_site_div_high[:, k] = (1/high_mut_rate) * mts.diversity([ts.samples()], windows=windows, mode='site')[:, 0]

# and for independent subsamples

sample_size = 100
assert(sample_size * nreps <= ts.num_samples)
sord = np.array(list(ts.samples()))
np.random.shuffle(sord)
sample_sets = []
for k in range(nreps):
    sample_sets.append(sord[k * sample_size + np.arange(sample_size)])

mts = msprime.mutate(ts, rate=mut_rate, keep=True)
sub_branch_div = ts.diversity(sample_sets, windows=windows, mode='branch')
sub_site_div = (1/mut_rate) * mts.diversity(sample_sets, windows=windows, mode='site')

maxy = 0.7 * max(np.max(np.abs(all_site_div)), np.max(np.abs(all_branch_div)),
                 np.max(np.abs(sub_site_div)), np.max(np.abs(sub_branch_div)))


######## Plot
fig = plt.figure(figsize=(6, 6))

# indep't mutations: low
ax = fig.add_subplot(311)
ax.set_ylim(0.0, maxy)
sl, = ax.plot(x, all_site_div[:,0],
             alpha=0.5,
             color='black',
             label="site")
sll = ax.plot(x, all_site_div[:,1:],
             alpha=0.5,
             color='black')
bl, = ax.plot(x, all_branch_div,
             alpha=1.0,
             color='red',
             label="branch")
vl = ax.axvline(spos[0], linestyle=":", label="sweep")
for k, (pos, sval) in enumerate(zip(spos, s)):
    if k > 0:
        vll = ax.axvline(pos, linestyle=":")
    # # text labels with selection coefficients?
    # at = ax.text(pos, 0.95*maxy, "s={:.2}".format(sval))

ax.set_xlabel("chromosome position (bp)")
ax.set_ylabel("mean TMRCA (π/μ)")
leg = ax.legend(
         loc = "upper left",
         bbox_to_anchor = (1.01, 1.01),
         frameon = False,
         borderpad=0)
atl = ax.text(ts.sequence_length, maxy, "μ={:.0}".format(low_mut_rate),
              horizontalalignment='right',
              verticalalignment='top')

# indep't mutations: high
ax = fig.add_subplot(312)
ax.set_ylim(0.0, maxy)
sl, = ax.plot(x, all_site_div_high[:,0],
             alpha=0.5,
             color='black',
             label="site")
sll = ax.plot(x, all_site_div_high[:,1:],
             alpha=0.5,
             color='black')
bl, = ax.plot(x, all_branch_div,
             alpha=1.0,
             color='red',
             label="branch")
vl = ax.axvline(spos[0], linestyle=":", label="sweep")
for k, (pos, sval) in enumerate(zip(spos, s)):
    if k > 0:
        vll = ax.axvline(pos, linestyle=":")
    # # text labels with selection coefficients?
    # at = ax.text(pos, 0.95*maxy, "s={:.2}".format(sval))

ax.set_xlabel("chromosome position (bp)")
ax.set_ylabel("mean TMRCA (π/μ)")
leg = ax.legend(
         loc = "upper left",
         bbox_to_anchor = (1.01, 1.01),
         frameon = False,
         borderpad=0)
atl = ax.text(ts.sequence_length, maxy, "μ={:.0}".format(high_mut_rate),
              horizontalalignment='right',
              verticalalignment='top')


# independant subsamples
ax = fig.add_subplot(313)
ax.set_ylim(0.0, maxy)
sl, = ax.plot(x, sub_site_div[:,0],
             alpha=0.5,
             color='black',
             label="site")
sll = ax.plot(x, sub_site_div[:,1:],
             alpha=0.5,
             color='black')
bl, = ax.plot(x, sub_branch_div[:,0],
             alpha=0.5,
             color='red',
             label="branch")
bll = ax.plot(x, sub_branch_div[:,1:],
             alpha=0.5,
             color='red')
vl = ax.axvline(spos[0], linestyle=":", label="sweep")
for k, (pos, sval) in enumerate(zip(spos, s)):
    if k > 0:
        vll = ax.axvline(pos, linestyle=":")
    # # text labels with selection coefficients?
    # at = ax.text(pos, 0.95*maxy, "s={:.2}".format(sval))

ax.set_xlabel("chromosome position (bp)")
ax.set_ylabel("mean TMRCA (π/μ)")
leg = ax.legend(
         loc = "upper left",
         bbox_to_anchor = (1.01, 1.01),
         frameon = False,
         borderpad=0)
atl = ax.text(ts.sequence_length, maxy, "μ={:.0}\n n={}".format(mut_rate, sample_size),
              horizontalalignment='right',
              verticalalignment='top')

# done
fig.savefig("{}.{}.diversity.pdf".format(outbase, mut_rate), bbox_inches = "tight")
plt.close(fig)


#### other plots

if False:
    # mutations at different rates
    mut_rates = [u * mut_rate for u in (0.1, 1, 10, 100)]
    mut_site_div = np.zeros((len(windows) - 1, len(mut_rates)))
    mts = ts
    for k in range(len(mut_rates)):
        rate = mut_rates[k]
        drate = rate - ([0.0] + mut_rates)[k]
        mts = msprime.mutate(ts, rate=drate, keep=True)
        mut_site_div[:, k] = (1/rate) * mts.diversity([ts.samples()], windows=windows, mode='site')[:, 0]

    # diff mut rates
    ax = fig.add_subplot(311)
    ax.set_ylim(0.0, maxy)
    mut_colors = plt.cm.brg(np.linspace(0, 0.4, len(mut_rates)))
    sll = [ax.plot(x, mut_site_div[:, k],
                   alpha=0.8,
                   color=mut_colors[k])[0] for k in range(len(mut_rates))]
    bl, = ax.plot(x, all_branch_div,
                 alpha=1.0,
                 color='red',
                 label="branch")
    vl = ax.axvline(spos[0], linestyle=":", label="sweep")
    for k, (pos, sval) in enumerate(zip(spos, s)):
        if k > 0:
            vll = ax.axvline(pos, linestyle=":")
        # # text labels with selection coefficients?
        # at = ax.text(pos, 0.95*maxy, "s={:.2}".format(sval))

    ax.set_xlabel("chromosome position (bp)")
    ax.set_ylabel("mean diversity (π)")
    leg = ax.legend(
             sll + [bl, vl],
             ["mu = {:.0}".format(u) for u in mut_rates] + ["branch", "sweep"],
             loc = "upper left",
             bbox_to_anchor = (1.01, 1.01),
             frameon = False,
             borderpad = 0)
