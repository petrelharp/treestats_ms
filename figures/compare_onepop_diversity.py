#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys


usage = """
Usage:
    {} (branch diversity tsv) (site diversity tsv)
Here `maskfile` should be a fasta of 0s and 1s.
""".format(sys.argv[0])

if len(sys.argv) != 3:
    raise ValueError(usage)


divfiles = { 'branch' : sys.argv[1],
             'site' : sys.argv[2] }

outbase = ".".join(divfiles['branch'].split(".")[:-2]) + f".site.ratio"

plotfile = f"{outbase}.pdf"

bdivs = { a : np.loadtxt(divfiles[a], skiprows=1) for a in divfiles }
num_windows = bdivs['branch'].shape[0]
windows = np.linspace(0, 63025522, num_windows + 1).astype('int')



bdivs['ratio'] = bdivs['site'] / bdivs['branch']

x = windows[:-1] + np.diff(windows)/2
x /= 10**6 # convert to Mb

fig = plt.figure(figsize=(7 * len(windows) / 64, 6))

for j, mode in enumerate(['site', 'branch', 'ratio']):
    ax = fig.add_subplot(3,1,j+1)
    ax.xaxis.set_tick_params(labelbottom=mode=="ratio")
    ax.set_title(mode.capitalize())
    ax.grid(axis="x")
    line, = ax.plot(x, bdivs[mode])

ax.set_xlabel("chr20 position (Megabases)")
fig.subplots_adjust(hspace=0.375)

fig.savefig(plotfile, bbox_inches = "tight")
plt.close(fig)
