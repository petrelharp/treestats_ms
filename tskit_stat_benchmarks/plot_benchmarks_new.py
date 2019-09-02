import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

yfactors=[0.4,0.3]
for infile,yf in zip(["benchmarks_including_copy.txt", "benchmarks_without_copy.txt"], yfactors):
    print(infile, yf)
    b = pd.read_csv(infile, sep=' ')
    b = b.round({'Nu':4})
    b['seconds_snp'] = b.seconds/b.nmutations
    grp = b.groupby(['toolkit'])
    fig = plt.figure()
    for n,g in grp:
        plt.loglog(g.nsam, g.seconds_snp,'o-',label=n)
    # p = sns.scatterplot(x='nsam', y='seconds_snp',
    #p = sns.scatterplot(x='nsam', y='seconds_snp',
    #                    hue='toolkit', style='Nu',
    #                    data=b, alpha= 0.5)
    # plot = p.get_figure()
    # #plt.ylim(0.0, 0.0006)
    # plt.xscale('log')
    # ylim=plt.gca().get_ylim()
    # ylim=(ylim[0]*1e-1,ylim[1]*yf)
    # plt.ylim(*ylim)
    plt.xlabel("Sample size (haploid)")
    plt.ylabel("Run time per variant (seconds)")
    #plt.legend(loc=2,bbox_to_anchor=(1,1))
    of = infile.replace("txt","pdf")
    plt.legend(loc='upper left',fontsize='x-small')
    plt.tight_layout()
    plt.savefig(of)
    plt.clf()
