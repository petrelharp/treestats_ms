import pandas as pd
import matplotlib.pyplot as plt

for infile in ["benchmarks_including_copy.txt", "benchmarks_without_copy.txt"]:
    b = pd.read_csv(infile, sep=' ')
    b = b.round({'Nu': 4})
    b['seconds_snp'] = b.seconds/b.nmutations
    grp = b.groupby(['toolkit'])
    fig = plt.figure()
    for n, g in grp:
        plt.loglog(g.nsam, g.seconds_snp, 'o-', label=n)
    plt.xlabel("Sample size (haploid)")
    plt.ylabel("Run time per variant (seconds)")
    of = infile.replace("txt", "pdf")
    plt.legend(loc='upper left', fontsize='x-small')
    plt.tight_layout()
    plt.savefig(of)
    plt.clf()
