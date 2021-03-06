#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys

# Roughly match the standard TGP color pallettes
tgp_region_colors = {
    "EAS": sns.color_palette("Greens", 2)[1],
    "EUR": sns.color_palette("Blues", 1)[0],
    "AFR": sns.color_palette("Wistia", 3)[0],
    "AMR": sns.color_palette("Reds", 2)[1],
    "SAS": sns.color_palette("Purples", 2)[1],
}

divfiles = { 'branch' : "1kg/relate_chr20.branch.1000000.diversity.tsv",
        'site' : "1kg/relate_chr20.site.1000000.diversity.tsv" }

plotfile = "1kg/relate_chr20_site_div_branch.1000000.diversity.pdf"

bdivs = { a : np.loadtxt(divfiles[a], skiprows=1) for a in divfiles }
num_windows = bdivs['branch'].shape[0]
windows = np.linspace(0, 63025522, num_windows + 1).astype('int')


colors = plt.get_cmap("tab20").colors

pop_names = ['CHB', 'JPT', 'CHS', 'CDX', 'KHV', 'CEU', 'TSI', 'FIN', 'GBR', 'IBS', 'YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB', 'MXL', 'PUR', 'CLM', 'PEL', 'GIH', 'PJL', 'BEB', 'STU', 'ITU']
superpops = ['EAS', 'EAS', 'EAS', 'EAS', 'EAS', 'EUR', 'EUR', 'EUR', 'EUR', 'EUR', 'AFR', 'AFR', 'AFR', 'AFR', 'AFR', 'AFR', 'AFR', 'AMR', 'AMR', 'AMR', 'AMR', 'SAS', 'SAS', 'SAS', 'SAS', 'SAS']
num_pops = len(pop_names)


supernames = ['EAS', 'EUR', 'AFR', 'AMR', 'SAS']
superpopcolors = tgp_region_colors
sbdivs = { a : np.zeros((bdivs[a].shape[0], len(supernames))) for a in bdivs }
for a in bdivs:
    nsp = [0 for _ in supernames]
    for j, p in enumerate(pop_names):
        k = supernames.index(superpops[j])
        sbdivs[a][:, k] = sbdivs[a][:, k] * k / (k + 1) + bdivs[a][:, j] / (k + 1)
        nsp[k] += 1

sbdivs['ratio'] = sbdivs['site'] / sbdivs['branch']

x = windows[:-1] + np.diff(windows)/2
x /= 10**6 # convert to Mb

fig = plt.figure(figsize=(7 * len(windows) / 64, 6))

lines = {}
for j, mode in enumerate(['site', 'branch', 'ratio']):
    ax = fig.add_subplot(3,1,j+1)
    ax.xaxis.set_tick_params(labelbottom=mode=="ratio")
    ax.set_title(mode.capitalize())
    ax.grid(axis="x")
    for j in range(sbdivs[mode].shape[1]):
        label = supernames[j]
        line, = ax.plot(
                x, sbdivs[mode][:, j],
                color=superpopcolors[supernames[j]])
        lines[label] = line

ax.set_xlabel("chr20 position (Megabases)")
print(lines)
names = sorted(lines.keys())
fig.legend(
    [lines[name] for name in names], names, loc="upper center", ncol=5)
fig.subplots_adjust(hspace=0.375)


fig.savefig(plotfile, bbox_inches = "tight")
plt.close(fig)
